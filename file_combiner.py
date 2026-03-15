#!/usr/bin/env python3
"""
Project File Combiner - Combines all project files into a single output file with a GUI
"""

import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import tkinter.ttk as ttk
from pathlib import Path
import os
from threading import Thread
import json


class FileCombinerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Project File Combiner")
        self.root.geometry("900x700")
        self.root.resizable(True, True)
        
        self.project_path = None
        self.blacklist = set()
        self.exclude_patterns = {
            '.git', '.gitignore', '__pycache__', '.pyc', '.pyo',
            'node_modules', '.env', '.venv', 'venv', 'dist', 'build',
            '.DS_Store', 'Thumbs.db', '.vscode', '.idea', '.settings'
        }
        
        self.setup_ui()
        self.load_config()
    
    def setup_ui(self):
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(4, weight=1)
        
        # Project Path Selection
        ttk.Label(main_frame, text="Project Path:", font=("Arial", 10, "bold")).grid(
            row=0, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        
        path_frame = ttk.Frame(main_frame)
        path_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))
        path_frame.columnconfigure(1, weight=1)
        
        self.path_var = tk.StringVar()
        path_entry = ttk.Entry(path_frame, textvariable=self.path_var, state='readonly')
        path_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))
        
        browse_btn = ttk.Button(path_frame, text="Browse...", command=self.select_project_path)
        browse_btn.grid(row=0, column=2)
        
        # Output Path Selection
        ttk.Label(main_frame, text="Output File:", font=("Arial", 10, "bold")).grid(
            row=2, column=0, columnspan=3, sticky=tk.W, pady=(0, 5)
        )
        
        output_frame = ttk.Frame(main_frame)
        output_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))
        output_frame.columnconfigure(1, weight=1)
        
        self.output_var = tk.StringVar(value="combined_output.txt")
        output_entry = ttk.Entry(output_frame, textvariable=self.output_var)
        output_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 5))
        
        output_browse_btn = ttk.Button(output_frame, text="Browse...", command=self.select_output_path)
        output_browse_btn.grid(row=0, column=2)
        
        # Blacklist Section
        ttk.Label(main_frame, text="Files/Folders to Exclude:", font=("Arial", 10, "bold")).grid(
            row=4, column=0, columnspan=3, sticky=tk.W, pady=(10, 5)
        )
        
        # Blacklist display area
        list_frame = ttk.Frame(main_frame)
        list_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(0, 10))
        list_frame.columnconfigure(0, weight=1)
        list_frame.rowconfigure(0, weight=1)
        
        scrollbar = ttk.Scrollbar(list_frame)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        self.blacklist_listbox = tk.Listbox(list_frame, height=10, yscrollcommand=scrollbar.set)
        self.blacklist_listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.config(command=self.blacklist_listbox.yview)
        
        # Blacklist buttons
        blacklist_btn_frame = ttk.Frame(main_frame)
        blacklist_btn_frame.grid(row=6, column=0, columnspan=3, sticky=tk.W, pady=(0, 10))
        
        ttk.Button(blacklist_btn_frame, text="Add to Blacklist", command=self.add_blacklist).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(blacklist_btn_frame, text="Remove Selected", command=self.remove_blacklist).pack(side=tk.LEFT)
        ttk.Button(blacklist_btn_frame, text="Clear All", command=self.clear_blacklist).pack(side=tk.LEFT, padx=(5, 0))
        
        # Options frame
        options_frame = ttk.LabelFrame(main_frame, text="Options", padding="10")
        options_frame.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))
        options_frame.columnconfigure(0, weight=1)
        
        self.include_binary_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(options_frame, text="Include binary files (images, etc.)", 
                       variable=self.include_binary_var).pack(anchor=tk.W)
        
        self.include_deps_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(options_frame, text="Include node_modules and dependencies", 
                       variable=self.include_deps_var).pack(anchor=tk.W)
        
        # Action buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=8, column=0, columnspan=3, sticky=(tk.W, tk.E))
        button_frame.columnconfigure(1, weight=1)
        
        ttk.Button(button_frame, text="Combine Files", command=self.combine_files).pack(side=tk.LEFT)
        ttk.Button(button_frame, text="Reset", command=self.reset_form).pack(side=tk.LEFT, padx=(5, 0))
        
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(button_frame, textvariable=self.status_var, foreground="blue").pack(side=tk.RIGHT)
    
    def select_project_path(self):
        path = filedialog.askdirectory(title="Select Project Root Directory")
        if path:
            self.project_path = Path(path)
            self.path_var.set(str(path))
            self.populate_files()
    
    def select_output_path(self):
        path = filedialog.asksaveasfilename(
            title="Save Combined File As",
            defaultextension=".txt",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        if path:
            self.output_var.set(path)
    
    def populate_files(self):
        """Populate the blacklist with default exclusions based on project contents"""
        if not self.project_path:
            return
        
        self.blacklist_listbox.delete(0, tk.END)
        
        # Add default exclusions
        for pattern in self.exclude_patterns:
            if (self.project_path / pattern).exists():
                self.blacklist_listbox.insert(tk.END, pattern)
    
    def add_blacklist(self):
        if not self.project_path:
            messagebox.showwarning("Warning", "Please select a project path first")
            return
        
        dialog = BlacklistDialog(self.root, self.project_path)
        self.root.wait_window(dialog.dialog)  # Wait for dialog to close
        if dialog.result:
            self.blacklist_listbox.insert(tk.END, dialog.result)
    
    def remove_blacklist(self):
        selection = self.blacklist_listbox.curselection()
        if selection:
            self.blacklist_listbox.delete(selection[0])
    
    def clear_blacklist(self):
        if messagebox.askyesno("Confirm", "Clear all blacklist entries?"):
            self.blacklist_listbox.delete(0, tk.END)
    
    def get_blacklist(self):
        """Get current blacklist from listbox"""
        items = self.blacklist_listbox.get(0, tk.END)
        return set(items)
    
    def should_exclude(self, file_path, blacklist):
        """Check if a file should be excluded"""
        relative_path = file_path.relative_to(self.project_path)
        
        # Check if any blacklist item matches
        for exclude_item in blacklist:
            # Check exact match or parent directory match
            if exclude_item in str(relative_path):
                return True
            # Check if it's a direct child
            if file_path.name == exclude_item:
                return True
        
        # Check file extension for binary files
        if not self.include_binary_var.get():
            binary_extensions = {
                '.png', '.jpg', '.jpeg', '.gif', '.bmp', '.ico', '.svg',
                '.mp3', '.mp4', '.wav', '.zip', '.rar', '.7z',
                '.exe', '.dll', '.so', '.o', '.pyc', '.pyo', '.class',
            }
            if file_path.suffix.lower() in binary_extensions:
                return True
        
        return False
    
    def combine_files(self):
        if not self.project_path:
            messagebox.showerror("Error", "Please select a project path")
            return
        
        output_path = self.output_var.get()
        if not output_path:
            messagebox.showerror("Error", "Please specify an output file")
            return
        
        self.status_var.set("Processing...")
        self.root.update()
        
        # Run in thread to keep UI responsive
        thread = Thread(target=self._do_combine, args=(output_path,))
        thread.daemon = True
        thread.start()
    
    def _do_combine(self, output_path):
        try:
            blacklist = self.get_blacklist()
            file_count = 0
            total_size = 0
            
            with open(output_path, 'w', encoding='utf-8', errors='replace') as output_file:
                # Walk through all files
                for file_path in sorted(self.project_path.rglob('*')):
                    if file_path.is_file():
                        # Skip blacklisted files
                        if self.should_exclude(file_path, blacklist):
                            continue
                        
                        try:
                            relative_path = file_path.relative_to(self.project_path)
                            
                            # Write file header
                            output_file.write(f"\n{'='*80}\n")
                            output_file.write(f"FILE: {relative_path}\n")
                            output_file.write(f"{'='*80}\n\n")
                            
                            # Try to read as text
                            try:
                                with open(file_path, 'r', encoding='utf-8') as f:
                                    content = f.read()
                                    output_file.write(content)
                            except UnicodeDecodeError:
                                output_file.write("[Binary file - not included]\n")
                                continue
                            
                            output_file.write("\n")
                            file_count += 1
                            total_size += file_path.stat().st_size
                        
                        except Exception as e:
                            output_file.write(f"[Error reading file: {e}]\n\n")
            
            # Success message
            self.status_var.set(f"✓ Combined {file_count} files ({total_size/1024:.1f} KB)")
            messagebox.showinfo(
                "Success",
                f"Files combined successfully!\n\n"
                f"Files included: {file_count}\n"
                f"Total size: {total_size/1024:.1f} KB\n"
                f"Output: {output_path}"
            )
        
        except Exception as e:
            self.status_var.set("Error occurred")
            messagebox.showerror("Error", f"An error occurred:\n{str(e)}")
    
    def reset_form(self):
        self.path_var.set("")
        self.output_var.set("combined_output.txt")
        self.project_path = None
        self.blacklist_listbox.delete(0, tk.END)
        self.status_var.set("Ready")
    
    def save_config(self):
        """Save user settings"""
        config = {
            'project_path': str(self.project_path) if self.project_path else None,
            'output_path': self.output_var.get(),
            'blacklist': list(self.blacklist_listbox.get(0, tk.END)),
            'include_binary': self.include_binary_var.get(),
            'include_deps': self.include_deps_var.get(),
        }
        try:
            with open('file_combiner_config.json', 'w') as f:
                json.dump(config, f, indent=2)
        except Exception:
            pass
    
    def load_config(self):
        """Load previous settings"""
        try:
            if os.path.exists('file_combiner_config.json'):
                with open('file_combiner_config.json', 'r') as f:
                    config = json.load(f)
                
                if config.get('project_path'):
                    self.project_path = Path(config['project_path'])
                    self.path_var.set(config['project_path'])
                    self.populate_files()
                
                self.output_var.set(config.get('output_path', 'combined_output.txt'))
                self.include_binary_var.set(config.get('include_binary', False))
                self.include_deps_var.set(config.get('include_deps', False))
                
                for item in config.get('blacklist', []):
                    self.blacklist_listbox.insert(tk.END, item)
        except Exception:
            pass


class BlacklistDialog:
    def __init__(self, parent, project_path):
        self.result = None
        self.project_path = project_path
        
        self.dialog = tk.Toplevel(parent)
        self.dialog.title("Add to Blacklist")
        self.dialog.geometry("400x300")
        self.dialog.grab_set()
        
        # Title
        ttk.Label(self.dialog, text="Select files/folders to exclude:", font=("Arial", 10, "bold")).pack(
            anchor=tk.W, padx=10, pady=(10, 5)
        )
        
        # File tree view
        frame = ttk.Frame(self.dialog)
        frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        
        scrollbar = ttk.Scrollbar(frame)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
        self.listbox = tk.Listbox(frame, yscrollcommand=scrollbar.set, height=12)
        self.listbox.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.config(command=self.listbox.yview)
        
        # Populate with top-level items
        try:
            for item in sorted(project_path.iterdir()):
                display_name = item.name
                if item.is_dir():
                    display_name += "/"
                self.listbox.insert(tk.END, display_name)
        except Exception:
            pass
        
        # Manual entry
        ttk.Label(self.dialog, text="Or enter custom pattern:", font=("Arial", 9)).pack(anchor=tk.W, padx=10, pady=(5, 0))
        
        entry_frame = ttk.Frame(self.dialog)
        entry_frame.pack(fill=tk.X, padx=10, pady=5)
        
        self.entry_var = tk.StringVar()
        ttk.Entry(entry_frame, textvariable=self.entry_var).pack(fill=tk.X)
        
        # Buttons
        button_frame = ttk.Frame(self.dialog)
        button_frame.pack(fill=tk.X, padx=10, pady=10)
        
        ttk.Button(button_frame, text="Add Selected", command=lambda: self.done(True)
                  ).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(button_frame, text="Add Custom", command=lambda: self.done(False)
                  ).pack(side=tk.LEFT)
        ttk.Button(button_frame, text="Cancel", command=lambda: self.done(None)
                  ).pack(side=tk.RIGHT)
    
    def done(self, use_selected):
        if use_selected is True:
            selection = self.listbox.curselection()
            if selection:
                self.result = self.listbox.get(selection[0]).rstrip('/')
        elif use_selected is False:
            custom = self.entry_var.get().strip()
            if custom:
                self.result = custom
        
        self.dialog.destroy()


def main():
    root = tk.Tk()
    app = FileCombinerApp(root)
    
    def on_closing():
        app.save_config()
        root.destroy()
    
    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()


if __name__ == "__main__":
    main()
