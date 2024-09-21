#!/usr/bin/env python3.12
# usage: python3.12 visual_codon.py

"""
This program calculates and plots the GC content of gene sequence,
so that the user can manually modify codons to achieve the goal of adjusting
local GC content, which is convenient for subsequent experiments,
such as subcloning or PCR gene synthesis.

The program needs Python3.12 and matplotlib 3.9.2.
"""

import sys
import os
from operator import itemgetter
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk


class App:
    def __init__(self, master):
        self.master = master

        # global variable
        self.gene_sequence = ''

        self.font_size_initial = 6
        self.x_label_spacing_initial = 20
        self.import_optimized_file = 0
        self.open_file_name = ""
        self.import_file_name = ""
        self.radiobutton_list = []
        self.selected_codon = tk.StringVar()
        # update codon flag
        self.updateflag = 0
        self.update_button = None
        self.label_aa_no_head = None
        self.label_aa_head = None
        self.label_original_codon_head = None
        self.label_optimized_codon_head = None
        self.label_blank_1 = None
        self.label_blank_2 = None
        self.label_aa_no = None
        self.label_aa = None
        self.label_original_codon = None
        self.label_optimized_codon = None
        self.label_blank_3 = None
        self.label_blank_4 = None
        self.radiobutton = None


        # dictionaries for amino acids and codons
        self.codon_to_amino_acid_dictionary = \
            {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TGT': 'C', 'TGC': 'C', 'GAT': 'D', 'GAC': 'D',
             'GAA': 'E', 'GAG': 'E', 'TTT': 'F', 'TTC': 'F', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
             'CAT': 'H', 'CAC': 'H', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'AAA': 'K', 'AAG': 'K', 'TTA': 'L',
             'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATG': 'M', 'AAT': 'N', 'AAC': 'N',
             'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
             'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
             'AGT': 'S', 'AGC': 'S', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GTT': 'V', 'GTC': 'V',
             'GTA': 'V', 'GTG': 'V', 'TGG': 'W', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGA': '*'
             }

        self.amino_acid_to_codons_with_frequency_dictionary = \
            {'A': [('GCG', 0.33), ('GCC', 0.26), ('GCA', 0.23), ('GCT', 0.18)],
             'C': [('TGC', 0.54), ('TGT', 0.46)],
             'D': [('GAT', 0.63), ('GAC', 0.37)],
             'E': [('GAA', 0.68), ('GAG', 0.32)],
             'F': [('TTT', 0.58), ('TTC', 0.42)],
             'G': [('GGC', 0.37), ('GGT', 0.35), ('GGG', 0.15), ('GGA', 0.13)],
             'H': [('CAT', 0.57), ('CAC', 0.43)],
             'I': [('ATT', 0.49), ('ATC', 0.39), ('ATA', 0.11)],
             'K': [('AAA', 0.74), ('AAG', 0.26)],
             'L': [('CTG', 0.47), ('TTA', 0.14), ('TTG', 0.13), ('CTT', 0.12), ('CTC', 0.1), ('CTA', 0.04)],
             'M': [('ATG', 1.0)],
             'N': [('AAC', 0.51), ('AAT', 0.49)],
             'P': [('CCG', 0.49), ('CCA', 0.2), ('CCT', 0.18), ('CCC', 0.13)],
             'Q': [('CAG', 0.66), ('CAA', 0.34)],
             'R': [('CGT', 0.36), ('CGC', 0.36), ('CGG', 0.11), ('CGA', 0.07), ('AGA', 0.07), ('AGG', 0.04)],
             'S': [('AGC', 0.25), ('TCT', 0.17), ('AGT', 0.16), ('TCC', 0.15), ('TCA', 0.14), ('TCG', 0.14)],
             'T': [('ACC', 0.4), ('ACG', 0.25), ('ACT', 0.19), ('ACA', 0.17)],
             'V': [('GTG', 0.35), ('GTT', 0.28), ('GTC', 0.2), ('GTA', 0.17)],
             'W': [('TGG', 1.0)],
             'Y': [('TAT', 0.59), ('TAC', 0.41)],
             '*': [('TAA', 0.61), ('TGA', 0.3), ('TAG', 0.09)]
             }

        # sort the codon frequencies for each amino acid
        for key in self.amino_acid_to_codons_with_frequency_dictionary.keys():
            self.amino_acid_to_codons_with_frequency_dictionary[key] = \
                    sorted(self.amino_acid_to_codons_with_frequency_dictionary[key], key=itemgetter(1), reverse=True)

        # amino_acid_to_codons_rank_dictionary
        self.amino_acid_to_codons_rank_dictionary = {}
        tmp_list = []
        for key in self.amino_acid_to_codons_with_frequency_dictionary.keys():
            for t in self.amino_acid_to_codons_with_frequency_dictionary[key]:
                tmp_list.append(t[0])
            self.amino_acid_to_codons_rank_dictionary[key] = tmp_list
            tmp_list = []

        master.title("Visual Codon V1.1")

        # get screen size
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        # root.geometry('1440x900+0+0') # for MacBook Air 13.3" with Retina Display
        master.geometry(str(screen_width) + 'x' + str(screen_height) + '+0+0')
        # master.wm_state('zoomed')
        master.protocol("WM_DELETE_WINDOW", self.on_closing)

        # create menu
        menu_bar = tk.Menu(root)
        master.config(menu=menu_bar)

        # Menu - File
        self.file_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="File", menu=self.file_menu)
        self.file_menu.add_command(label="Open Gene 1 to edit", command=self.open_file)
        self.file_menu.add_command(label="Save optimized gene", command=self.save_optimized_gene)
        self.file_menu.add_command(label="Export table to txt", command=self.export_table_to_txt)
        self.file_menu.add_command(label="Export changed codons", command=self.export_changed_codons_to_txt)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Import Gene 2 to compare", command=self.import_second_gene_sequence)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Reinitialize", command=self.initialize_program)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command=self.exit_program)

        # Menu - Edit
        edit_menu = tk.Menu(menu_bar, tearoff=0)
        edit_menu.add_command(label="Set graph X axis", command=self.set_graph_X_axis)
        menu_bar.add_cascade(label="Edit", menu=edit_menu)

        # create frame1
        frame1 = ttk.Frame(root, width=900, height=300)
        frame1.pack(side="top", fill="both", expand=True)

        # create TreeView with scrollbar
        self.tree = ttk.Treeview(frame1, columns=("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8"),
                                 selectmode="browse", show="headings")

        scrollbar = ttk.Scrollbar(frame1, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)

        # create frame2
        self.frame2 = ttk.Frame(root, width=900, height=100)
        self.frame2.pack(side="top", pady=5)

        # add TreeView headers
        self.tree.heading("col1", text="amino acid No")
        self.tree.heading("col2", text="amino acid")
        self.tree.heading("col3", text="original codon")
        self.tree.heading("col4", text="original rank")
        self.tree.heading("col5", text="original GC%")
        self.tree.heading("col6", text="optimized codon")
        self.tree.heading("col7", text="optimized rank")
        self.tree.heading("col8", text="optimized GC%")

        # column width
        self.tree.column("col1", width=100, anchor="center")
        self.tree.column("col2", width=100, anchor="center")
        self.tree.column("col3", width=100, anchor="center")
        self.tree.column("col4", width=100, anchor="center")
        self.tree.column("col5", width=100, anchor="center")
        self.tree.column("col6", width=100, anchor="center")
        self.tree.column("col7", width=100, anchor="center")
        self.tree.column("col8", width=100, anchor="center")

        # create frame3
        frame3 = ttk.Frame(root, width=900, height=500)

        frame3.pack(side="top", fill="both", pady=5, expand=True)

        # create output area for Matplotlib graph
        self.figure = plt.figure(figsize=(6, 2))
        self.canvas = FigureCanvasTkAgg(self.figure, master=frame3)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)

        toolbar = NavigationToolbar2Tk(self.canvas, frame3)
        toolbar.update()
        toolbar.pack(side="top", fill=tk.X)

        self.ax = self.figure.add_subplot(111)
        plt.axis('off')

    @staticmethod
    def read_sequence_from_file(file_path):
        """
        read sequence file
        """
        sequences = []
        current_sequence = ''
        with open(file_path, 'r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                if not line.startswith('>'):
                    current_sequence += line
                else:
                    if current_sequence:
                        sequences.append(current_sequence)
                        current_sequence = ''
            if current_sequence:
                sequences.append(current_sequence)
            sequence_str = ''.join([i.upper() for i in sequences if i.isalpha()])
        return sequence_str
    def open_file(self):
        """
        function to open and check gene sequence
        """
        file_path = filedialog.askopenfilename()
        if file_path:
            if file_path[-6:] == '.fasta' or file_path[-3:] == '.fa':
                self.initialize_program()
                self.gene_sequence = ''
                self.import_optimized_file = 0
                self.updateflag = 0
                self.open_file_name = os.path.basename(file_path)
                self.gene_sequence = self.read_sequence_from_file(file_path)
                print("The original gene sequence is as follows.")
                print("\n", end='')
                print(self.gene_sequence)
                print("\n", end='')

                print("The length of the gene sequence is %s." % len(self.gene_sequence))
                print("\n", end='')

                # test gene length
                if len(self.gene_sequence) < 60:
                    print("The length of gene sequence is less than 60.")
                    print("Please synthesize directly.")
                    messagebox.showerror(
                        "Error",
                        f"The length of gene sequence is less than 60. Please synthesize directly."
                    )
                    return

                # test codons
                if len(self.gene_sequence) % 3 != 0:
                    print("The length of gene sequence is not the folds of three.")
                    print("Please check the gene fasta file. Thank you!")
                    messagebox.showerror(
                        "Error",
                        f"The length of gene sequence is not the folds of three. Please check the gene fasta file."
                    )
                    return

                # test if there is any stop codon before the stop codon in the end
                for i in range(int(len(self.gene_sequence)/3)-1):
                    codon = self.gene_sequence[3*i:3*i+3]
                    if codon in {'TAA', 'TAG', 'TGA'}:
                        position = str(3*i + 1) + ',' + str(3*i + 2) + ',' + str(3*i + 3)
                        print("There is a stop codon in the sequence before the stop codon.")
                        print("The position is " + position + ".")
                        print("Please check the gene fasta file. Thank you!")
                        messagebox.showerror(
                            "Error",
                            f"There is a stop codon in the sequence before the stop codon. The position is {position}. "
                            f"Please check the gene fasta file."
                        )
                        return

                # test if all bases are ATCG
                for i in range(len(self.gene_sequence)):
                    if self.gene_sequence[i] not in {'A', 'T', 'C', 'G'}:
                        print("There is a base other than ATCG.")
                        print("It is %s in position %s." % (str(self.gene_sequence[i]), str(i+1)))
                        print("Please check the gene fasta file. Thank you!")
                        messagebox.showerror(
                            "Error",
                            f"There is a base other than ATCG. It is {str(self.gene_sequence[i])} in position {str(i+1)}. "
                            f"Please check the gene fasta file."
                        )
                        return

                # fill in the table for the first time
                self.insert_items()

                # left mouse single click
                self.tree.bind('<ButtonRelease-1>', self.select_row)

                # change window title
                self.master.title("Visual Codon V1.1: " + file_path)
                self.file_menu.entryconfig(5, state="normal")
                self.update_button.config(state=tk.NORMAL)
            else:
                print("It is not a fasta file.")
                print("Please re-select a fasta file.")
                messagebox.showerror(
                    "Error",
                    f"It is not a fasta file. Please re-select a fasta file."
                )
                return

    def save_optimized_gene(self):
        """
        function to save the optimized gene sequence to a fasta file
        """
        if self.gene_sequence:
            file_path = filedialog.asksaveasfilename()
            if file_path:
                if file_path[-6:] == '.fasta':
                    file_path = file_path[0:-6]
                file = open(file_path + '.fasta', 'w')
                file.write('>optimized gene sequence\n')
                file.write(self.gene_sequence + '\n')
                file.close()

    def export_table_to_txt(self):
        """
        function to export table to a txt file.
        The txt file can be opened and copied to macOS Numbers. The data are all in cells.
        """
        if self.gene_sequence:
            file_path = filedialog.asksaveasfilename()
            if file_path:
                if file_path[-4:] == '.txt':
                    file_path = file_path[0:-4]
                file = open(file_path + '.txt', 'w')
                # combine the words of each header
                header = [self.tree.heading(col)["text"].replace(" ", "_") for col in self.tree["columns"]]
                # write headers
                file.write("\t".join(header) + "\n")
                # write codon rows
                for i, item in enumerate(self.tree.get_children()):
                    values = [str(self.tree.item(item)["text"])]
                    for col in self.tree["columns"]:
                        if col == "col4" or col == "col7":
                            values.append("\'" + str(self.tree.item(item)["values"][self.tree["columns"].index(col)]))
                        else:
                            values.append(str(self.tree.item(item)["values"][self.tree["columns"].index(col)]))
                    # [1:] to delete the space on the first left position for each data rows
                    file.write("\t".join(values)[1:] + "\n")
                file.close()

    def export_changed_codons_to_txt(self):
        """
        function to export changed codons to a txt file.
        The txt file can be opened and copied to macOS Numbers. The data are all in cells.
        """
        if self.gene_sequence:
            file_path = filedialog.asksaveasfilename()
            if file_path:
                if file_path[-4:] == '.txt':
                    file_path = file_path[0:-4]
                file = open(file_path + '.txt', 'w')
                # combine the words of each header
                header = [self.tree.heading(col)["text"].replace(" ", "_") for col in self.tree["columns"]]
                # write headers
                file.write("\t".join(header) + "\n")
                # write codon rows
                for i, item in enumerate(self.tree.get_children()):
                    values = [str(self.tree.item(item)["text"])]
                    for col in self.tree["columns"]:
                        if self.tree.item(item)["values"][self.tree["columns"].index("col4")] != self.tree.item(item)["values"] \
                            [self.tree["columns"].index("col7")]:
                            if col == "col4" or col == "col7":
                                values.append("\'" + str(self.tree.item(item)["values"][self.tree["columns"].index(col)]))
                            else:
                                values.append(str(self.tree.item(item)["values"][self.tree["columns"].index(col)]))

                    # [1:] to delete the space on the first left position for each data rows
                    line = "\t".join(values)[1:] + "\n"
                    if not line.isspace():
                        file.write(line)
                file.close()

    def import_second_gene_sequence(self):
        """
        function to open another gene sequence to compare.
        """
        # if self.updateflag == 1:
        #     messagebox.showinfo("warning", "")
        imported_optimized_gene_sequence = ''
        # test if gene 1 is opened
        if len(self.gene_sequence) == 0:
            print("Please open gene 1 first. Thank you!")
            messagebox.showerror(
                "Error",
                f"Please open gene 1 first."
            )
            return

        file_path = filedialog.askopenfilename()
        # print("optimized file_path is " + file_path)
        if file_path:
            if file_path[-6:] == '.fasta' or file_path[-3:] == '.fa':
                self.import_file_name = os.path.basename(file_path)
                imported_optimized_gene_sequence = self.read_sequence_from_file(file_path)
                print("The imported optimized gene sequence is as follows.")
                print("\n", end='')
                print(imported_optimized_gene_sequence)
                print("\n", end='')

                print("The length of the imported optimized gene sequence is %s." % len(imported_optimized_gene_sequence))
                print("\n", end='')

                # test gene length
                if len(imported_optimized_gene_sequence) < 60:
                    print("The length of the imported optimized gene sequence is less than 60.")
                    print("Please synthesize directly.")
                    messagebox.showerror(
                        "Error",
                        f"The length of the imported optimized gene sequence is less than 60. Please synthesize directly. "
                    )
                    return

                # test codons
                if len(imported_optimized_gene_sequence) % 3 != 0:
                    print("The length of the imported optimized gene sequence is not the folds of three.")
                    print("Please check the gene fasta file. Thank you!")
                    messagebox.showerror(
                        "Error",
                        f"The length of the imported optimized gene sequence is not the folds of three. "
                        f"Please check the gene fasta file. "
                    )
                    return

                # test if there is any stop codon before the stop codon in the end
                for i in range(int(len(imported_optimized_gene_sequence)/3)-1):
                    codon = imported_optimized_gene_sequence[3*i:3*i+3]
                    if codon in {'TAA', 'TAG', 'TGA'}:
                        position = str(3*i + 1) + ',' + str(3*i + 2) + ',' + str(3*i + 3)
                        print("There is a stop codon in the sequence before the stop codon.")
                        print("The position is " + position + ".")
                        print("Please check the gene fasta file. Thank you!")
                        messagebox.showerror(
                            "Error",
                            f"There is a stop codon in the sequence before the stop codon. "
                            f"The position is {position}. Please check the gene fasta file."
                        )
                        return

                # test if all bases are ATCG
                for i in range(len(imported_optimized_gene_sequence)):
                    if imported_optimized_gene_sequence[i] not in {'A', 'T', 'C', 'G'}:
                        print("There is a base other than ATCG.")
                        print("It is %s in position %s." % (str(imported_optimized_gene_sequence[i]), str(i+1)))
                        print("Please check the gene fasta file. Thank you!")
                        messagebox.showerror(
                            "Error",
                            f"There is a base other than ATCG. It is {str(imported_optimized_gene_sequence[i])} "
                            f"in position {str(i+1)}. Please check the gene fasta file."
                        )
                        return

                # test if length of two genes are the same
                if len(imported_optimized_gene_sequence) != len(self.gene_sequence):
                    print("Length of two genes are different. Please check the gene fasta file. Thank you!")
                    messagebox.showerror(
                        "Error",
                        f"Length of two genes are different. Length of gene 1={len(self.gene_sequence)} "
                        f"Length of gene 2={len(imported_optimized_gene_sequence)}. Please check the imported gene fasta file."
                    )
                    return

                # compare amino acid sequences of the original gene sequence and the imported optimized gene sequence
                # if the original gene has been opened
                if self.gene_sequence:
                    amino_acid_sequence_of_original_gene = ''
                    amino_acid_sequence_of_optimized_gene = ''
                    for i in range(int(len(self.gene_sequence)/3)):
                        codon = self.gene_sequence[3*i:3*i+3]
                        amino_acid = self.codon_to_amino_acid_dictionary[codon]
                        amino_acid_sequence_of_original_gene = amino_acid_sequence_of_original_gene + amino_acid

                    for i in range(int(len(imported_optimized_gene_sequence)/3)):
                        codon = imported_optimized_gene_sequence[3*i:3*i+3]
                        amino_acid = self.codon_to_amino_acid_dictionary[codon]
                        amino_acid_sequence_of_optimized_gene = amino_acid_sequence_of_optimized_gene + amino_acid

                    if amino_acid_sequence_of_original_gene == amino_acid_sequence_of_optimized_gene:
                        # fill in the columns 6 and column 7
                        for i in range(int(len(imported_optimized_gene_sequence)/3)):
                            row_id = self.tree.get_children()[i]
                            codon = imported_optimized_gene_sequence[3*i:3*i+3]
                            codon_rank = str(self.amino_acid_to_codons_rank_dictionary[self.codon_to_amino_acid_dictionary[codon]].
                                             index(codon)+1) + '/' + str(len(self.amino_acid_to_codons_rank_dictionary
                                                                             [self.codon_to_amino_acid_dictionary[codon]]))
                            self.tree.set(row_id, "col6", codon)
                            self.tree.set(row_id, "col7", codon_rank)

                        # fill in the column 8
                        for i in range(3, int(len(imported_optimized_gene_sequence)/3)-3):
                            subsequence = ''
                            for j in range(i-3,i+4):
                                row_id = self.tree.get_children()[j]
                                subsequence = subsequence + self.tree.item(row_id, option="values")[5]
                            GC_percentage = round((subsequence.count('G') + subsequence.count('C')) * 100 / len(subsequence), 2)
                            row_id = self.tree.get_children()[i]
                            self.tree.set(row_id, "col8", GC_percentage)

                        # set the background color of the row if the original codon and the optimized codon are different
                        for i in range(0,int(len(imported_optimized_gene_sequence)/3)):
                            row_id = self.tree.get_children()[i]
                            if self.tree.item(row_id, option="values")[2] != self.tree.item(row_id, option="values")[5]:
                                # define style and apply to the row
                                self.tree.tag_configure("new_bg", background='LightGray')
                                self.tree.item(row_id, tags=("new_bg",))
                            else:
                                self.tree.tag_configure("white_bg", background='White')
                                self.tree.item(row_id, tags=("white_bg",))

                        self.import_optimized_file = 1
                        self.update_gc_graph(self.x_label_spacing_initial, self.font_size_initial)
                        print("The optimized gene sequence has been imported and the table updated.")
                        self.update_button.config(state=tk.DISABLED)
                        self.master.title("Visual Codon V1.1: comparison mode [ReadOnly]")
                else:
                    print("Please open the original gene fasta file first.")
            else:
                print("It is not a fasta file.")
                print("Please re-select a fasta file.")

    def initialize_program(self):
        """
        reinitialize this program
        :return:
        """
        self.master.title("Visual Codon V1.1")
        self.file_menu.entryconfig(5, state="normal")
        # global variable
        self.gene_sequence = ''

        self.font_size_initial = 6
        self.x_label_spacing_initial = 20
        self.import_optimized_file = 0
        self.open_file_name = ""
        self.import_file_name = ""
        # self.radiobutton_list.clear()
        # self.selected_codon = tk.StringVar()
        # update codon flag
        self.updateflag = 0
        self.clear_center_display()

        self.tree.delete(*self.tree.get_children())
        self.ax.cla()
        self.canvas.draw_idle()
        plt.axis('off')

    def clear_center_display(self):
        self.update_button and self.update_button.grid_forget()
        for radiobutton in self.radiobutton_list:
            radiobutton and radiobutton.grid_forget()
        self.label_aa_no_head and self.label_aa_no_head.grid_forget()
        self.label_aa_head and self.label_aa_head.grid_forget()
        self.label_original_codon_head and self.label_original_codon_head.grid_forget()
        self.label_optimized_codon_head and self.label_optimized_codon_head.grid_forget()
        self.label_blank_1 and self.label_blank_1.grid_forget()
        self.label_blank_2 and self.label_blank_2.grid_forget()
        self.label_aa_no and self.label_aa_no.grid_forget()
        self.label_aa and self.label_aa.grid_forget()
        self.label_original_codon and self.label_original_codon.grid_forget()
        self.label_optimized_codon and self.label_optimized_codon.grid_forget()
        self.label_blank_3 and self.label_blank_3.grid_forget()
        self.label_blank_4 and self.label_blank_4.grid_forget()
        # self.radiobutton and self.radiobutton.grid_forget()


    def exit_program(self):
        """
        function to exit.
        """
        root.destroy()
        sys.exit()

    def update_selected_item(self):
        """
        function to update the selected tree row.
        """
        selected_item = self.tree.selection()
        if selected_item:
            # change codon
            self.tree.set(selected_item, "col6", self.selected_codon.get())
            new_rank = str(self.amino_acid_to_codons_rank_dictionary[self.codon_to_amino_acid_dictionary[self.
                           selected_codon.get()]].index(self.selected_codon.get())+1) + '/' + \
                           str(len(self.amino_acid_to_codons_rank_dictionary[self.codon_to_amino_acid_dictionary[self.selected_codon.get()]]))

            self.tree.set(selected_item, "col7", new_rank)

            # change gene_sequence
            # gene_codon_number = int(len(self.gene_sequence)/3)
            # self.gene_sequence = ''
            # for i in range(gene_codon_number):
            #     row_id = self.tree.get_children()[i]
            #     self.gene_sequence = self.gene_sequence + self.tree.item(row_id, option="values")[5]

            # change optimized GC%
            row = self.tree.item(selected_item[0])['values'][0] - 1
            # only change 7 rows
            for i in range(row-3, row+4):
                if i in range(3, int(len(self.gene_sequence) / 3) - 3):
                    subsequence = ''
                    for j in range(i - 3, i + 4):
                        row_id = self.tree.get_children()[j]
                        subsequence = subsequence + self.tree.item(row_id, option="values")[5]
                    GC_percentage = round((subsequence.count('G') + subsequence.count('C')) * 100 / len(subsequence), 2)
                    row_id = self.tree.get_children()[i]
                    self.tree.set(row_id, "col8", GC_percentage)

            # set the background color of the row if codon changed
            if self.tree.item(self.tree.selection()[0])['values'][2] != self.tree.item(self.tree.selection()[0])['values'][5]:
                # define style and apply to the row
                self.tree.tag_configure("new_bg", background='LightGray')
                self.tree.item(self.tree.selection()[0], tags=("new_bg",))
            else:
                self.tree.tag_configure("white_bg", background='White')
                self.tree.item(self.tree.selection()[0], tags=("white_bg",))

            self.update_gc_graph(self.x_label_spacing_initial, self.font_size_initial)

    def select_row(self, event):
        """
        function for single left mouse click.
        """
        # get the selected row
        selected_item = self.tree.selection()
        if selected_item:
            item = selected_item[0]     # The selected_item is a tuple, add [0] to get the first element of the tuple.
            values = self.tree.item(item)['values']
            # print(f"Selected row: {values}")

        self.clear_center_display()
        # labels - headers
        self.label_aa_no_head = tk.Label(self.frame2, text="amino acid No", bg='LightGray', width=14, height=1)
        self.label_aa_no_head.grid(row=0, column=0, padx=1, pady=1)

        self.label_aa_head = tk.Label(self.frame2, text="amino acid", bg='LightGray', width=14, height=1)
        self.label_aa_head.grid(row=0, column=1, padx=1, pady=1)

        self.label_original_codon_head = tk.Label(self.frame2, text="original codon", bg='LightGray', width=14, height=1)
        self.label_original_codon_head.grid(row=0, column=2, padx=1, pady=1)

        self.label_optimized_codon_head = tk.Label(self.frame2, text="optimized codon", bg='LightGray', width=14, height=1)
        self.label_optimized_codon_head.grid(row=0, column=3, padx=1, pady=1)

        # blank labels - headers
        self.label_blank_1 = tk.Label(self.frame2, text="   ", width=14, height=1)
        self.label_blank_1.grid(row=0, column=4, padx=1, pady=1)

        self.label_blank_2 = tk.Label(self.frame2, text="   ", width=14, height=1)
        self.label_blank_2.grid(row=0, column=5, padx=1, pady=1)

        # labels - data
        self.label_aa_no = tk.Label(self.frame2, text=values[0], bg='White', width=14, height=1)
        self.label_aa_no.grid(row=1, column=0, padx=1, pady=1)

        self.label_aa = tk.Label(self.frame2, text=values[1], bg='White', width=14, height=1)
        self.label_aa.grid(row=1, column=1, padx=1, pady=1)

        self.label_original_codon = tk.Label(self.frame2, text=values[2], bg='White', width=14, height=1)
        self.label_original_codon.grid(row=1, column=2, padx=1, pady=1)

        self.label_optimized_codon = tk.Label(self.frame2, text=values[5], bg='White', width=14, height=1)
        self.label_optimized_codon.grid(row=1, column=3, padx=1, pady=1)

        # blank labels - data
        self.label_blank_3 = tk.Label(self.frame2, text="   ", width=14, height=1)
        self.label_blank_3.grid(row=1, column=4, padx=1, pady=1)

        self.label_blank_4 = tk.Label(self.frame2, text="   ", width=14, height=1)
        self.label_blank_4.grid(row=1, column=5, padx=1, pady=1)

        # clear the previous radiobutton
        # for radiobutton in self.radiobutton_list:
        #     radiobutton.destroy()

        # create radiobutton(s)
        Radio_button_number = len(self.amino_acid_to_codons_with_frequency_dictionary[values[1]])
        for i in range(Radio_button_number):
            eligible_codon_frequency = self.amino_acid_to_codons_with_frequency_dictionary[values[1]][i]

            # delete '
            self.radiobutton = tk.Radiobutton(self.frame2, text=f"{eligible_codon_frequency[0]}, {eligible_codon_frequency[1]}", \
                                         variable=self.selected_codon, value=eligible_codon_frequency[0])
            self.radiobutton.grid(row=2, column=i, padx=1, pady=1)
            self.radiobutton_list.append(self.radiobutton)

        # label - blank
        # for i in range(Radio_button_number, 6):
        #
        #     self.label_optimized_codon_head = tk.Label(self.frame2, text="   ", width=14, height=1)
        #     self.label_optimized_codon_head.grid(row=2, column=i, padx=1, pady=1)

        # default radiobutton
        self.selected_codon.set(values[2])

        # create update codon button
        self.update_button = ttk.Button(self.frame2, text="Update codon", command=self.update_selected_item)
        self.update_button.grid(row=2, column=6, padx=1, pady=1)
        if self.import_optimized_file == 1:
            self.update_button.config(state=tk.DISABLED)
            self.master.title("Visual Codon V1.1: comparison mode [ReadOnly]")
        else:
            self.update_button.config(state=tk.NORMAL)

    def set_X_axis_func(self, set_graph_X_axis, font_size_entry, font_spacing_entry):
        """
        function to set x axis value.
        command.
        """
        font_size_str = font_size_entry.get()
        font_spacing_str = font_spacing_entry.get()
        if font_size_str and font_spacing_str:
            try:
                font_size = int(font_size_str)
                x_label_spacing = int(font_spacing_str)
                self.update_gc_graph(x_label_spacing, font_size)
                set_graph_X_axis.destroy()
            except ValueError:
                messagebox.showerror("Error", "Please enter a valid integer!")
        else:
            messagebox.showerror("Error", "Please enter integers!")

    def set_graph_X_axis(self):
        """
        function to set x axis value.
        """
        if self.gene_sequence:

            set_graph_X_axis = tk.Toplevel(root)
            set_graph_X_axis.title("set_graph_X_axis")

            font_size_label = tk.Label(set_graph_X_axis, text="Please input font size for X axis:")
            font_size_label.pack(side="top")

            font_size_entry = tk.Entry(set_graph_X_axis)
            font_size_entry.insert(0, "6")  # default '5'
            font_size_entry.pack(side="top", padx=10, pady=5)

            font_spacing_label = tk.Label(set_graph_X_axis, text="Please input spacing for X axis:")
            font_spacing_label.pack(side="top")

            font_spacing_entry = tk.Entry(set_graph_X_axis)
            font_spacing_entry.insert(0, "20")  # default '20'
            font_spacing_entry.pack(side="top", padx=10, pady=5)

            ok_button = tk.Button(set_graph_X_axis, text="OK", command=lambda: self.set_X_axis_func(set_graph_X_axis, font_size_entry, font_spacing_entry))
            ok_button.pack(side="top", pady=5)


    def on_mouse_press(self, event):
        """
        function to get right row in treeview.
        """
        if event.inaxes:
            x, y = event.xdata, event.ydata
            # coord_label.config(text=f"Mouse Coordinates: ({x:.2f}, {y:.2f})")
            intx = int(x)
            if intx in range(0, int(len(self.gene_sequence) / 3)):
                self.tree.selection_set(self.tree.get_children()[intx])
                self.tree.yview(self.tree.index(self.tree.get_children()[intx]))
                self.select_row(None)

    def insert_items(self):
        """
        function to fill in the table when open the fasta file
        need to delete the rows, insert new rows and update the graph in frame3
        """
        # delete the existing rows
        self.tree.delete(*self.tree.get_children())
        # fill in columns 1-4
        for i in range(int(len(self.gene_sequence)/3)):
            codon = self.gene_sequence[3*i:3*i+3]
            amino_acid = self.codon_to_amino_acid_dictionary[codon]
            rank = str(self.amino_acid_to_codons_rank_dictionary[amino_acid].index(codon)+1) + \
                   '/' + str(len(self.amino_acid_to_codons_rank_dictionary[amino_acid]))

            self.tree.insert("", "end", values=(i+1, amino_acid, codon, rank))

        # fill in the colomn 5
        for i in range(3, int(len(self.gene_sequence)/3)-3):
            subsequence = ''
            for j in range(i-3,i+4):
                row_id = self.tree.get_children()[j]
                subsequence = subsequence + self.tree.item(row_id, option="values")[2]

            GC_percentage = round((subsequence.count('G') + subsequence.count('C')) * 100 / len(subsequence), 2)

            row_id = self.tree.get_children()[i]
            self.tree.set(row_id, "col5", GC_percentage)

        # fill in the columns 6-7
        for i in range(int(len(self.gene_sequence)/3)):
            row_id = self.tree.get_children()[i]
            original_codon = self.tree.item(row_id, option="values")[2]
            original_rank = self.tree.item(row_id, option="values")[3]

            self.tree.set(row_id, "col6", original_codon)
            self.tree.set(row_id, "col7", original_rank)

        # fill in the column 8
        for i in range(3, int(len(self.gene_sequence)/3)-3):
            row_id = self.tree.get_children()[i]
            original_gc_percentage = self.tree.item(row_id, option="values")[4]
            self.tree.set(row_id, "col8", original_gc_percentage)

        # set default selected row
        self.tree.selection_set(self.tree.get_children()[0])

        # click the 1st row (by default). The previous line has selected the first row. Now do something, i.e., select_row.
        self.select_row(None)

        # auto adjust x label spacing

        self.x_label_spacing_initial = int(int(len(self.gene_sequence)/3) / 20)
        # print(x_label_spacing_initial)
        self.update_gc_graph(self.x_label_spacing_initial, self.font_size_initial)

    def on_closing(self):
        """
        function to destroy canvas before exit program
        """
        self.canvas.get_tk_widget().destroy()
        self.exit_program()

    def update_gc_graph(self, x_label_spacing_in, font_size_in):
        """
        function to update the Matplotlib gc graph
        """
        plt.axis('on')
        codon_no_list = [i for i in range(3+1, int(len(self.gene_sequence)/3)-3+1)]
        original_gc_percentage_list = []
        optimized_gc_percentage_list = []

        for i in range(3,int(len(self.gene_sequence)/3)-3):
            row_id = self.tree.get_children()[i]
            original_gc_percentage = self.tree.item(row_id, option="values")[4]
            optimized_gc_percentage = self.tree.item(row_id, option="values")[7]

            original_gc_percentage_list.append(original_gc_percentage)
            optimized_gc_percentage_list.append(optimized_gc_percentage)

        self.ax.clear()

        plt.xlabel('Codon No.', fontsize=6)
        plt.ylabel('GC percentage', fontsize=6)

        # set the distance between the axis label and the coordinate axis
        self.ax.xaxis.set_label_coords(0.5, -0.1)
        self.ax.yaxis.set_label_coords(-0.05, 0.5)

        # X axis label from 1 to len
        # print(f"x_label_spacing_in={x_label_spacing_in}")
        plt.xticks(range(1, int(len(self.gene_sequence)/3)+1, x_label_spacing_in), fontsize=font_size_in)
        plt.yticks(range(int(min(original_gc_percentage_list)), int(max(original_gc_percentage_list))+1, 5), fontsize=font_size_in)

        if self.import_optimized_file == 0:
            self.ax.plot(codon_no_list, original_gc_percentage_list, linewidth=0.3, label="Gene 1:"+self.open_file_name, color='brown')
            if optimized_gc_percentage_list != original_gc_percentage_list:
                self.ax.plot(codon_no_list, optimized_gc_percentage_list, linewidth=0.3, label="Gene 1 optimized", color='b')
                self.updateflag = 1
                self.file_menu.entryconfig(5, state=tk.DISABLED)
        else:
            self.ax.plot(codon_no_list, original_gc_percentage_list, linewidth=0.3, label="Gene 1:"+self.open_file_name, color='brown')
            self.ax.plot(codon_no_list, optimized_gc_percentage_list, linewidth=0.3, label="Gene 2:"+self.import_file_name, color='b')
        self.ax.legend(loc='upper right', frameon=False, fontsize=6)

        # draw line of GC percentages of original gene sequence
        sequence = ''
        for i in range(int(len(self.gene_sequence)/3)):
            row_id = self.tree.get_children()[i]
            sequence = sequence + self.tree.item(row_id, option="values")[2]

        GC_percentage = round((sequence.count('G') + sequence.count('C')) * 100 / len(sequence), 2)
        self.ax.axhline(y=GC_percentage, c='brown', ls='--', lw=0.2)

        # draw line of GC percentages of optimized gene sequence
        if self.import_optimized_file == 1 or optimized_gc_percentage_list != original_gc_percentage_list:
            sequence = ''
            for i in range(int(len(self.gene_sequence)/3)):
                row_id = self.tree.get_children()[i]
                sequence = sequence + self.tree.item(row_id, option="values")[5]

            GC_percentage = round((sequence.count('G') + sequence.count('C')) * 100 / len(sequence), 2)
            self.ax.axhline(y=GC_percentage, c='b', ls='--', lw=0.2)

        self.canvas.draw()

        cid = self.canvas.mpl_connect('button_press_event', self.on_mouse_press)


if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()

