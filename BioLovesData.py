#!/usr/bin/env python3
"""
BioLovesData - Gene Analysis Program
Integrated Bioinformatics Analysis Tool using Tkinter GUI
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
import sys

# Import analysis modules
from calculate_transversion_ratio import calculate_transversion_ratio_analysis
from consensus_symbols import consensus_symbols_analysis
from find_promoter import find_promoter_analysis
from find_tm_domain import find_tm_domain_analysis
from multi_seq_alignment import multi_seq_alignment_analysis

class BioLovesDataGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("BioLovesData - Gene Analysis Program")
        self.root.geometry("1200x800")
        self.root.configure(bg='#f0f0f0')
        
        # Currently selected file
        self.current_file = None
        self.current_file_content = ""
        
        # Initialize GUI components
        self.setup_gui()
        
    def setup_gui(self):
        """Setup GUI components"""
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Title
        title_label = ttk.Label(main_frame, text="BioLovesData", font=("Arial", 20, "bold"))
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        # Left panel (File selection and function selection)
        left_frame = ttk.LabelFrame(main_frame, text="Analysis Settings", padding="10")
        left_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 10))
        
        # File selection section
        file_frame = ttk.LabelFrame(left_frame, text="File Selection", padding="5")
        file_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        
        self.file_label = ttk.Label(file_frame, text="Selected File: None", wraplength=250)
        self.file_label.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 5))
        
        ttk.Button(file_frame, text="Select File", command=self.select_file).grid(row=1, column=0, sticky=(tk.W, tk.E), padx=(0, 5))
        ttk.Button(file_frame, text="Preview File", command=self.preview_file).grid(row=1, column=1, sticky=(tk.W, tk.E))
        
        # Analysis function selection section
        function_frame = ttk.LabelFrame(left_frame, text="Analysis Function", padding="5")
        function_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        
        self.analysis_var = tk.StringVar(value="transversion")
        
        analyses = [
            ("Transversion Ratio", "transversion", "Calculate transition/transversion"),
            ("Consensus Symbols", "consensus", "Analyze conservation of MSA"),
            ("Promoter Finder", "promoter", "Detect promoter regions"),
            ("TM Domain Finder", "tm_domain", "Predict transmembrane domains"),
            ("Multi-Seq Alignment", "alignment", "Perform MSA")
        ]
        
        for i, (name, value, desc) in enumerate(analyses):
            rb = ttk.Radiobutton(function_frame, text=name, variable=self.analysis_var, value=value)
            rb.grid(row=i, column=0, sticky=tk.W, pady=2)
            
            desc_label = ttk.Label(function_frame, text=desc, font=("Arial", 8), foreground="gray")
            desc_label.grid(row=i, column=1, sticky=tk.W, padx=(10, 0))
        
        # Options setting section
        self.options_frame = ttk.LabelFrame(left_frame, text="Analysis Options", padding="5")
        self.options_frame.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        
        # Run button
        ttk.Button(left_frame, text="Run Analysis", command=self.run_analysis, 
                  style="Accent.TButton").grid(row=3, column=0, sticky=(tk.W, tk.E), pady=10)
        
        # Right panel (Result display)
        right_frame = ttk.LabelFrame(main_frame, text="Analysis Results", padding="10")
        right_frame.grid(row=1, column=1, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Result text area
        self.result_text = scrolledtext.ScrolledText(right_frame, width=70, height=35, wrap=tk.WORD)
        self.result_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Save results button
        ttk.Button(right_frame, text="Save Results", command=self.save_results).grid(row=1, column=0, pady=(10, 0))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(1, weight=1)
        left_frame.columnconfigure(0, weight=1)
        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)
        
        # Update options when analysis function changes
        self.analysis_var.trace('w', self.update_options)
        self.update_options()
        
    def select_file(self):
        """File selection dialog"""
        file_types = [
            ("All Supported Files", "*.fasta;*.fa;*.txt"),
            ("FASTA Files", "*.fasta;*.fa"),
            ("Text Files", "*.txt"),
            ("All Files", "*.*")
        ]
        
        filename = filedialog.askopenfilename(
            title="Select a file for analysis",
            filetypes=file_types
        )
        
        if filename:
            self.current_file = filename
            try:
                with open(filename, 'r', encoding='utf-8') as f:
                    self.current_file_content = f.read()
                
                # Display only filename
                display_name = os.path.basename(filename)
                if len(display_name) > 30:
                    display_name = display_name[:27] + "..."
                
                self.file_label.config(text=f"Selected File: {display_name}")
                
            except Exception as e:
                messagebox.showerror("Error", f"Could not read file: {str(e)}")
                self.current_file = None
                self.current_file_content = ""
    
    def preview_file(self):
        """File preview"""
        if not self.current_file_content:
            messagebox.showwarning("Warning", "Please select a file first.")
            return
        
        preview_window = tk.Toplevel(self.root)
        preview_window.title("File Preview")
        preview_window.geometry("600x400")
        
        preview_text = scrolledtext.ScrolledText(preview_window, wrap=tk.WORD)
        preview_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Display only first 1000 characters
        content = self.current_file_content[:1000]
        if len(self.current_file_content) > 1000:
            content += "\n\n... (File is longer)"
        
        preview_text.insert(tk.END, content)
        preview_text.config(state=tk.DISABLED)
    
    def update_options(self, *args):
        """Update options based on selected analysis"""
        # Remove existing option widgets
        for widget in self.options_frame.winfo_children():
            widget.destroy()
        
        analysis_type = self.analysis_var.get()
        
        if analysis_type == "transversion":
            ttk.Label(self.options_frame, text="No additional options.").grid(row=0, column=0)
            
        elif analysis_type == "consensus":
            ttk.Label(self.options_frame, text="Input Format:").grid(row=0, column=0, sticky=tk.W)
            self.consensus_format = tk.StringVar(value="fasta")
            ttk.Radiobutton(self.options_frame, text="FASTA File", 
                           variable=self.consensus_format, value="fasta").grid(row=1, column=0, sticky=tk.W)
            ttk.Radiobutton(self.options_frame, text="Space-separated Sequences", 
                           variable=self.consensus_format, value="space").grid(row=2, column=0, sticky=tk.W)
            
        elif analysis_type == "promoter":
            ttk.Label(self.options_frame, text="Sequence Length Limit:").grid(row=0, column=0, sticky=tk.W)
            self.promoter_length = tk.StringVar(value="")
            ttk.Entry(self.options_frame, textvariable=self.promoter_length, width=10).grid(row=0, column=1, sticky=tk.W)
            ttk.Label(self.options_frame, text="(Leave empty for full sequence)").grid(row=0, column=2, sticky=tk.W)
            
            ttk.Label(self.options_frame, text="Max Ambiguity:").grid(row=1, column=0, sticky=tk.W)
            self.promoter_ambiguity = tk.StringVar(value="0.001")
            ttk.Entry(self.options_frame, textvariable=self.promoter_ambiguity, width=10).grid(row=1, column=1, sticky=tk.W)
            
        elif analysis_type == "tm_domain":
            ttk.Label(self.options_frame, text="Hydrophobicity Threshold:").grid(row=0, column=0, sticky=tk.W)
            self.tm_threshold = tk.StringVar(value="1.5")
            ttk.Entry(self.options_frame, textvariable=self.tm_threshold, width=10).grid(row=0, column=1, sticky=tk.W)
            
            ttk.Label(self.options_frame, text="Hydrophobicity Scale:").grid(row=1, column=0, sticky=tk.W)
            self.tm_scale = tk.StringVar(value="KyteDoolittle")
            scale_combo = ttk.Combobox(self.options_frame, textvariable=self.tm_scale, width=15)
            scale_combo['values'] = ["KyteDoolittle", "WimleyWhite", "HessaHigy", "MoonFleming", "TanfordTaylor"]
            scale_combo.grid(row=1, column=1, sticky=tk.W)
            scale_combo.state(['readonly'])
            
        elif analysis_type == "alignment":
            ttk.Label(self.options_frame, text="Scoring Matrix:").grid(row=0, column=0, sticky=tk.W)
            self.scoring_matrix = tk.StringVar(value="BLOSUM62")
            matrix_combo = ttk.Combobox(self.options_frame, textvariable=self.scoring_matrix, width=20)
            matrix_combo['values'] = ["BLOSUM62", "TEAM_SAMWISE", "defaultScoringMatrix", "HoxD55", "HoxD70"]
            matrix_combo.grid(row=0, column=1, sticky=tk.W)
            matrix_combo.state(['readonly'])
            
            ttk.Label(self.options_frame, text="Gap Penalty:").grid(row=1, column=0, sticky=tk.W)
            self.gap_penalty = tk.StringVar(value="-2")
            ttk.Entry(self.options_frame, textvariable=self.gap_penalty, width=10).grid(row=1, column=1, sticky=tk.W)
    
    def run_analysis(self):
        """Run selected analysis"""
        if not self.current_file_content:
            messagebox.showwarning("Warning", "Please select a file first.")
            return
        
        analysis_type = self.analysis_var.get()
        
        try:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, "Analyzing...\n")
            self.root.update()
            
            if analysis_type == "transversion":
                result = calculate_transversion_ratio_analysis(self.current_file_content)
                self.display_transversion_result(result)
                
            elif analysis_type == "consensus":
                if hasattr(self, 'consensus_format') and self.consensus_format.get() == "space":
                    # Process as space-separated sequences
                    sequences = self.current_file_content.strip().split()
                    result = consensus_symbols_analysis(sequences)
                else:
                    # Process as FASTA file
                    from calculate_transversion_ratio import read_seq_file_content
                    info_arr, seq_arr = read_seq_file_content(self.current_file_content)
                    result = consensus_symbols_analysis(seq_arr)
                self.display_consensus_result(result)
                
            elif analysis_type == "promoter":
                seq_length = None
                if hasattr(self, 'promoter_length') and self.promoter_length.get():
                    try:
                        seq_length = int(self.promoter_length.get())
                    except ValueError:
                        pass
                
                max_ambiguity = 0.001
                if hasattr(self, 'promoter_ambiguity') and self.promoter_ambiguity.get():
                    try:
                        max_ambiguity = float(self.promoter_ambiguity.get())
                    except ValueError:
                        pass
                
                result = find_promoter_analysis(self.current_file_content, seq_length, max_ambiguity)
                self.display_promoter_result(result)
                
            elif analysis_type == "tm_domain":
                threshold = 1.5
                if hasattr(self, 'tm_threshold') and self.tm_threshold.get():
                    try:
                        threshold = float(self.tm_threshold.get())
                    except ValueError:
                        pass
                
                scale = "KyteDoolittle"
                if hasattr(self, 'tm_scale') and self.tm_scale.get():
                    scale = self.tm_scale.get()
                
                result = find_tm_domain_analysis(self.current_file_content, threshold, scale)
                self.display_tm_domain_result(result)
                
            elif analysis_type == "alignment":
                gap_penalty = -2
                if hasattr(self, 'gap_penalty') and self.gap_penalty.get():
                    try:
                        gap_penalty = int(self.gap_penalty.get())
                    except ValueError:
                        pass
                
                scoring_matrix_name = "BLOSUM62"
                if hasattr(self, 'scoring_matrix') and self.scoring_matrix.get():
                    scoring_matrix_name = self.scoring_matrix.get()
                
                result = multi_seq_alignment_analysis(self.current_file_content, gap_penalty=gap_penalty, scoring_matrix_name=scoring_matrix_name)
                self.display_alignment_result(result)
                
        except Exception as e:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, f"An error occurred during analysis:\n{str(e)}")
    
    def display_transversion_result(self, result):
        """Display Transversion ratio result"""
        self.result_text.delete(1.0, tk.END)
        
        if "error" in result:
            self.result_text.insert(tk.END, f"Error: {result['error']}\n")
            return
        
        self.result_text.insert(tk.END, "=== Transversion Ratio Analysis Result ===\n\n")
        self.result_text.insert(tk.END, f"Number of sequences analyzed: {result['sequences_count']}\n")
        self.result_text.insert(tk.END, f"Transversions: {result['transversions']}\n")
        self.result_text.insert(tk.END, f"Transitions: {result['transitions']}\n")
        self.result_text.insert(tk.END, f"Transversion/Transition Ratio: {result['ratio']}\n\n")
        self.result_text.insert(tk.END, f"Result: {result['message']}\n")
    
    def display_consensus_result(self, result):
        """Display Consensus symbols result"""
        self.result_text.delete(1.0, tk.END)
        
        if "error" in result:
            self.result_text.insert(tk.END, f"Error: {result['error']}\n")
            return
        
        self.result_text.insert(tk.END, "=== Consensus Symbols Analysis Result ===\n\n")
        
        # Display sequences
        for i, seq in enumerate(result['sequences']):
            self.result_text.insert(tk.END, f"Sequence {i+1}: {seq}\n")
        
        self.result_text.insert(tk.END, f"Symbols: {result['symbols']}\n\n")
        
        # Statistics
        self.result_text.insert(tk.END, "=== Statistics ===\n")
        self.result_text.insert(tk.END, f"Total sites: {result['total_sites']}\n")
        for symbol, count in result['symbol_counts'].items():
            self.result_text.insert(tk.END, f"{symbol}: {count}\n")
        
        self.result_text.insert(tk.END, "\n=== Legend ===\n")
        for symbol, meaning in result['legend'].items():
            self.result_text.insert(tk.END, f"'{symbol}': {meaning}\n")
    
    def display_promoter_result(self, result):
        """Display Promoter finder result"""
        self.result_text.delete(1.0, tk.END)
        
        if "error" in result:
            self.result_text.insert(tk.END, f"Error: {result['error']}\n")
            return
        
        self.result_text.insert(tk.END, "=== Promoter Finder Analysis Result ===\n\n")
        self.result_text.insert(tk.END, f"Sequence Description: {result['description']}\n")
        self.result_text.insert(tk.END, f"Analyzed Sequence Length: {result['sequence_length']}\n")
        self.result_text.insert(tk.END, f"Number of Promoters Found: {result['total_promoters']}\n\n")
        
        if result['total_promoters'] > 0:
            for i, promoter in enumerate(result['promoter_regions']):
                self.result_text.insert(tk.END, f"--- Promoter {i+1} ---\n")
                self.result_text.insert(tk.END, f"-10 box location: {promoter['box_10_location']}\n")
                self.result_text.insert(tk.END, f"-10 box sequence: {promoter['box_10_sequence']} (Similarity: {promoter['box_10_similarity']})\n")
                self.result_text.insert(tk.END, f"Spacer length: {promoter['spacer_length']}\n")
                self.result_text.insert(tk.END, f"Spacer sequence: {promoter['spacer_sequence']}\n")
                self.result_text.insert(tk.END, f"-35 box location: {promoter['box_35_location']}\n")
                self.result_text.insert(tk.END, f"-35 box sequence: {promoter['box_35_sequence']} (Similarity: {promoter['box_35_similarity']})\n")
                self.result_text.insert(tk.END, f"Full promoter: {promoter['full_promoter']}\n\n")
        
        # Parameter information
        params = result['parameters']
        self.result_text.insert(tk.END, "=== Analysis Parameters ===\n")
        self.result_text.insert(tk.END, f"-10 box consensus: {params['box_10']}\n")
        self.result_text.insert(tk.END, f"-35 box consensus: {params['box_35']}\n")
        self.result_text.insert(tk.END, f"Similarity threshold: {params['threshold']}\n")
        self.result_text.insert(tk.END, f"Spacer length range: {params['min_spacer_len']}-{params['max_spacer_len']}\n")
    
    def display_tm_domain_result(self, result):
        """Display TM domain finder result"""
        self.result_text.delete(1.0, tk.END)
        
        if "error" in result:
            self.result_text.insert(tk.END, f"Error: {result['error']}\n")
            return
        
        self.result_text.insert(tk.END, "=== TM Domain Finder Analysis Result ===\n\n")
        
        for seq_result in result['results']:
            self.result_text.insert(tk.END, f"Sequence: {seq_result['name']}\n")
            
            if not seq_result['valid']:
                self.result_text.insert(tk.END, f"Error: {seq_result['error']}\n\n")
                continue
            
            self.result_text.insert(tk.END, f"Number of TM domains: {seq_result['tm_domain_count']}\n\n")
            
            if seq_result['tm_domain_count'] > 0:
                for i, domain in enumerate(seq_result['tm_domains']):
                    self.result_text.insert(tk.END, f"--- TM Domain {i+1} ---\n")
                    self.result_text.insert(tk.END, f"Location: {domain['location']}\n")
                    self.result_text.insert(tk.END, f"Length: {domain['length']}\n")
                    self.result_text.insert(tk.END, f"Hydrophobicity Value: {domain['hydrophobicity']}\n")
                    self.result_text.insert(tk.END, f"Sequence: {domain['sequence']}\n\n")
                
                # Display annotated sequence
                self.result_text.insert(tk.END, "=== Annotated Sequence ===\n")
                for segment in seq_result['annotated_sequence']:
                    if segment['type'] == 'tm_domain':
                        self.result_text.insert(tk.END, f"[TM:{segment['hydrophobicity']}]{segment['sequence']}")
                    else:
                        self.result_text.insert(tk.END, segment['sequence'])
                self.result_text.insert(tk.END, "\n\n")
            
        # Parameter information
        params = result['parameters']
        self.result_text.insert(tk.END, "=== Analysis Parameters ===\n")
        self.result_text.insert(tk.END, f"Hydrophobicity Threshold: {params['threshold']}\n")
        self.result_text.insert(tk.END, f"Hydrophobicity Scale: {params['scale']}\n")
        self.result_text.insert(tk.END, f"Window Size: {params['min_window']}-{params['max_window']}\n")
    
    def display_alignment_result(self, result):
        """Display Multiple sequence alignment result"""
        self.result_text.delete(1.0, tk.END)
        
        if "error" in result:
            self.result_text.insert(tk.END, f"Error: {result['error']}\n")
            return
        
        self.result_text.insert(tk.END, "=== Multiple Sequence Alignment Result ===\n\n")
        self.result_text.insert(tk.END, f"Number of Sequences: {result['sequence_count']}\n")
        self.result_text.insert(tk.END, f"Alignment Length: {result['alignment_length']}\n\n")
        
        # Display aligned sequences
        self.result_text.insert(tk.END, "=== Aligned Sequences ===\n")
        for seq_info in result['sequences']:
            self.result_text.insert(tk.END, f"{seq_info['name'][:20]:20} {seq_info['aligned_sequence']}\n")
        
        # Gap statistics
        gap_stats = result['gap_statistics']
        self.result_text.insert(tk.END, f"\n=== Gap Statistics ===\n")
        self.result_text.insert(tk.END, f"Number of positions with gaps: {gap_stats['total_gap_positions']}\n")
        
        if gap_stats['gap_positions']:
            self.result_text.insert(tk.END, "\nKey Gap Positions (first 20):\n")
            for gap_pos in gap_stats['gap_positions']:
                self.result_text.insert(tk.END, f"Position {gap_pos['position']}: {gap_pos['gap_count']} sequences ({gap_pos['gap_percentage']}%)\n")
        
        # Parameter information
        params = result['parameters']
        self.result_text.insert(tk.END, f"\n=== Analysis Parameters ===\n")
        self.result_text.insert(tk.END, f"Scoring Matrix: {params['scoring_matrix']}\n")
        self.result_text.insert(tk.END, f"Gap Penalty: {params['gap_penalty']}\n")
    
    def save_results(self):
        """Save results to file"""
        content = self.result_text.get(1.0, tk.END)
        if not content.strip():
            messagebox.showwarning("Warning", "No results to save.")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Save Results",
            defaultextension=".txt",
            filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
        )
        
        if filename:
            try:
                with open(filename, 'w', encoding='utf-8') as f:
                    f.write(content)
                messagebox.showinfo("Success", f"Results saved to: {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred while saving the file: {str(e)}")

def main():
    """Main function"""
    root = tk.Tk()
    
    # Style settings
    style = ttk.Style()
    style.theme_use('clam')
    
    app = BioLovesDataGUI(root)
    
    # Handle program exit
    def on_closing():
        root.quit()
        root.destroy()
    
    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()

if __name__ == "__main__":
    main()

