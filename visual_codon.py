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
import re

class InitialDialog(tk.Toplevel):
    def __init__(self, master):
        super().__init__(master)
        self.title("Select a host organism")
        self.geometry("400x250")

        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        window_width = 400
        window_height = 250

        x = int((screen_width / 2) - (window_width / 2))
        y = int((screen_height / 2) - (window_height / 2))

        self.geometry(f"{window_width}x{window_height}+{x}+{y}")

        self.label = tk.Label(self, text="Please select a host organism below:")
        self.label.pack(pady=5)

        # from https://www.genscript.com/tools/codon-frequency-table
        self.combobox = ttk.Combobox(self, values=["Escherichia coil", "Yeast", "Insect", "C. elegans", \
                                                   "Drosophila melanogaster", "Human", "Mouse", "Rat", \
                                                   "Pig", "Pichia pastoris", "Arabidopsis thaliana", \
                                                   "Streptomyces", "Zea mays (Maize)", \
                                                   "Nicotiana tabacum (Tabacco)", \
                                                   "Saccharomyces cerevisiae (gbpln)", \
                                                   "Cricetulus griseus (CHO)"])
        self.combobox.pack(pady=5)
        self.combobox.current(0)

        self.combobox.bind("<<ComboboxSelected>>", self.on_combobox_selected)

        self.selected_option = self.combobox.get()
        self.selected_filepath = None

        self.file_path = tk.StringVar()
        self.file_path.set("No file selected")

        self.label1 = tk.Label(self, text="OR load a host organism codon table file:")
        self.label1.pack(pady=5)

        self.open_button = tk.Button(self, text="Open File...", command=self.open_file)
        self.open_button.pack(pady=5)

        self.label2 = tk.Label(self, text="", justify='left', anchor='w', foreground='red')
        self.label2.pack(padx=3, pady=5)

        self.confirm_button = tk.Button(self, text="Confirm", command=self.confirm)
        self.confirm_button.pack(pady=5)

    def open_file(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.file_path.set(file_path)
            self.selected_filepath = file_path
            new_filepath = '\n'.join([file_path[i:i + 60] for i in range(0, len(file_path), 60)])
            self.label2.config(text=new_filepath)

    def on_combobox_selected(self, event):
        self.selected_option = self.combobox.get()

    def confirm(self):
        selected_file = self.file_path.get()
        if selected_file != "No file selected":
            print(f"Selected file: {selected_file}")
            # can do other check by parsing lines
            lines = []
            with open(self.selected_filepath, 'r', encoding='utf-8') as file:
                line_count = 0
                for line in file:
                    line = line.strip()
                    if line:
                        line_count += 1
                        lines.append(line)

                # can do other check, including line_count
                if line_count != 32:
                    messagebox.showerror(
                        "Error",
                        f"The line count of {self.selected_filepath} is {line_count} which should be 32."  
                        f" Please check the file."
                    )
                    return

            self.destroy()
        else:
            print(f"Selected file: None")
            self.destroy()
            # messagebox.showwarning("Warning", "Please select a file first!")

    def destroy(self):
        if self.selected_option is not None or self.selected_filepath is not None:
            # send message to master
            print("destroy initial dialog")
            self.master.on_dialog_closed(self.selected_option, self.selected_filepath)
        super().destroy()

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.codon_to_amino_acid_dictionary = {}
        self.amino_acid_to_codons_with_frequency_dictionary = {}
        self.stop_codons = []

        self.withdraw()  # hide
        self.dialog = InitialDialog(self)
        self.wait_window(self.dialog)  # wait
        self.deiconify()  # show

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

        self.title("Visual Codon V1.2")

        # get screen size
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        # root.geometry('1440x900+0+0') # for MacBook Air 13.3" with Retina Display
        self.geometry(str(screen_width) + 'x' + str(screen_height) + '+0+0')
        # master.wm_state('zoomed')
        self.protocol("WM_DELETE_WINDOW", self.on_closing)

        # create menu
        menu_bar = tk.Menu(self)
        self.config(menu=menu_bar)

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
        edit_menu.add_command(label="Adjust graph font", command=self.set_graph_X_axis)
        edit_menu.add_command(label="Customize organism", command=self.customize_organism)
        menu_bar.add_cascade(label="Edit", menu=edit_menu)

        # create frame1
        frame1 = ttk.Frame(self, width=900, height=300)
        frame1.pack(side="top", fill="both", expand=True)

        # create TreeView with scrollbar
        self.tree = ttk.Treeview(frame1, columns=("col1", "col2", "col3", "col4", "col5", "col6", "col7", "col8"),
                                 selectmode="browse", show="headings")

        scrollbar = ttk.Scrollbar(frame1, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        self.tree.pack(side="left", fill="both", expand=True)

        # create frame2
        self.frame2 = ttk.Frame(self, width=900, height=100)
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
        frame3 = ttk.Frame(self, width=900, height=500)

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

    def on_dialog_closed(self, selected_option, selected_filepath):
        # This function is called when the InitialDialog is closed and it receives two parameters.
        self.selected_option = selected_option
        self.selected_filepath = selected_filepath
        print(self.selected_option)
        print(self.selected_filepath)

        data_str = ""
        lines = []

        if self.selected_filepath is None:
            # from https://www.genscript.com/tools/codon-frequency-table
            # "Escherichia coil", "Yeast", "Insect", "C. elegans", "Drosophila melanogaster", \
            #     "Human", "Mouse", "Rat", "Pig", "Pichia pastoris", "Arabidopsis thaliana", \
            #     "Streptomyces", "Zea mays (Maize)", "Nicotiana tabacum (Tabacco)", "Saccharomyces cerevisiae (gbpln)", \
            #     "Cricetulus griseus (CHO)"
            if self.selected_option == "Escherichia coil":
                data_str = App.Ec
            elif self.selected_option == "Yeast":
                data_str = App.Yeast
            elif self.selected_option == "Insect":
                data_str = App.Insect
            elif self.selected_option == "C. elegans":
                data_str = App.Cele
            elif self.selected_option == "Drosophila melanogaster":
                data_str = App.Droso
            elif self.selected_option == "Human":
                data_str = App.Human
            elif self.selected_option == "Mouse":
                data_str = App.Mouse
            elif self.selected_option == "Rat":
                data_str = App.Rat
            elif self.selected_option == "Pig":
                data_str = App.Pig
            elif self.selected_option == "Pichia pastoris":
                data_str = App.Pichia
            elif self.selected_option == "Arabidopsis thaliana":
                data_str = App.Arabid
            elif self.selected_option == "Streptomyces":
                data_str = App.Strep
            elif self.selected_option == "Zea mays (Maize)":
                data_str = App.Zea
            elif self.selected_option == "Nicotiana tabacum (Tabacco)":
                data_str = App.Nicotia
            elif self.selected_option == "Saccharomyces cerevisiae (gbpln)":
                data_str = App.Saccha
            elif self.selected_option == "Cricetulus griseus (CHO)":
                data_str = App.Cricet

            lines = [line.strip() for line in data_str.split('\n') if line.strip()]
        else:
            with open(self.selected_filepath, 'r', encoding='utf-8') as file:
                line_count = 0
                for line in file:
                    line = line.strip()
                    if line:
                        line_count += 1
                        lines.append(line)

                # can do other check, besides line_count
                if line_count != 32:
                    messagebox.showerror(
                        "Error",
                        f"The line count of {self.selected_filepath} is {line_count} which should be 32. Please check the file."
                    )
                    print("Use Escherichia coil as default.")
                    self.selected_option = "Escherichia coil"
                    self.selected_filepath = ""
                    lines = []
                    data_str = ""
                    data_str = App.Ec
                    lines = [line.strip() for line in data_str.split('\n') if line.strip()]
                else:
                    self.selected_option = self.selected_filepath

        print(f"lines is {lines}")

        # parse lines
        codon_dict = {}
        codon_frequency_dict = {}
        stop_codon_list = []

        for line in lines:
            parts = line.split()
            codon_dict[parts[0]] = parts[1]
            codon_dict[parts[3]] = parts[4]

            # generate acid dictionary
            if parts[1] not in codon_frequency_dict:
                codon_frequency_dict[parts[1]] = [(parts[0], parts[2])]
            else:
                codon_frequency_dict[parts[1]].append((parts[0], parts[2]))
            if parts[4] not in codon_frequency_dict:
                codon_frequency_dict[parts[4]] = [(parts[3], parts[5])]
            else:
                codon_frequency_dict[parts[4]].append((parts[3], parts[5]))

            if parts[1] == '*':
                stop_codon_list.append(parts[0])
            if parts[4] == '*':
                stop_codon_list.append(parts[3])

        print(codon_dict)
        # print(codon_frequency_dict)
        for amino_acid, codon_list in codon_frequency_dict.items():
            print(f"{amino_acid}:{codon_list}")

        print(f"stop codons: {stop_codon_list}")
        self.codon_to_amino_acid_dictionary = codon_dict
        self.amino_acid_to_codons_with_frequency_dictionary = codon_frequency_dict
        self.stop_codons = stop_codon_list

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

                print(f"The length of the gene sequence is {len(self.gene_sequence)}.")
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
                    if codon in self.stop_codons:
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
                        print(f"It is {str(self.gene_sequence[i])} in position {str(i+1)}.")
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
                self.title("Visual Codon V1.2: " + file_path + " (\"" + self.selected_option + "\")")
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

                print(f"The length of the imported optimized gene sequence is {len(imported_optimized_gene_sequence)}.")
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
                    if codon in self.stop_codons:
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
                        print(f"It is {str(imported_optimized_gene_sequence[i])} in position {str(i+1)}.")
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
                        self.title("Visual Codon V1.2: comparison mode [ReadOnly]" + " (\"" + self.selected_option + "\")")
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
        info = self.selected_filepath if self.selected_filepath else self.selected_option
        messagebox.showinfo("info", f"The host organism is {info}. If you need to change the host organism, please restart this program.")


        self.title("Visual Codon V1.2")
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
        self.destroy()
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
            self.title("Visual Codon V1.2: comparison mode [ReadOnly]")
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

            set_graph_X_axis = tk.Toplevel(self)
            set_graph_X_axis.title("Adjust graph font")

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

    def customize_organism(self):
        """
        function to customize organism.
        """
        # definition of codons table
        codons = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC',
                  'TGA', 'TGG',
                  'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC',
                  'CGA', 'CGG',
                  'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC',
                  'AGA', 'AGG',
                  'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC',
                  'GGA', 'GGG']
        # codons = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC',
        #           'TGA', 'TGG']
        # codons = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG']
        customize_organism = tk.Toplevel(self)
        customize_organism.title("Customize a host organism codon table")
        screen_width = customize_organism.winfo_screenwidth()
        screen_height = customize_organism.winfo_screenheight()
        print(f"customize_organism screen_width={screen_width}, screen_height={screen_height}")
        # create table frame
        table_frame = tk.Frame(customize_organism)
        table_frame.pack()

        # head
        header_labels = ['Codon', 'Amino acid', 'Fraction', 'Codon', 'Amino acid', 'Fraction']
        for i, header_text in enumerate(header_labels):
            header_label = tk.Label(table_frame, text=header_text, font=("Arial", 8), relief=tk.FLAT, width=15)
            header_label.grid(row=0, column=i)

        # table content
        colors = ["lightgray", "white"]
        for row in range(len(codons)//2):
            color = colors[(row // 4) % 2]
            bg_color = color

            # codon1_label = tk.Label(table_frame, text=codons[8*(row // 4)+row % 4], font=("Arial", 10), relief=tk.RIDGE, width=15, bg=bg_color)
            codon1_label = tk.Label(table_frame, text=codons[8 * (row // 4) + row % 4], font=("Arial", 8), relief=tk.FLAT,
                                    width=15, bg=bg_color)
            codon1_label.grid(row=row + 1, column=0, pady=(0, 0))
            amino_acid1_entry = tk.Entry(table_frame, font=("Arial", 8), relief=tk.FLAT, width=15)
            amino_acid1_entry.grid(row=row + 1, column=1, pady=(0, 0))
            fraction1_entry = tk.Entry(table_frame, font=("Arial", 8), relief=tk.FLAT, width=15)
            fraction1_entry.grid(row=row + 1, column=2, pady=(0, 0))
            # codon2_label = tk.Label(table_frame, text=codons[8*(row // 4)+row % 4 + 4], font=("Arial", 10), relief=tk.RIDGE, width=15, bg=bg_color)
            codon2_label = tk.Label(table_frame, text=codons[8 * (row // 4) + row % 4 + 4], font=("Arial", 8), relief=tk.FLAT,
                                    width=15, bg=bg_color)
            codon2_label.grid(row=row + 1, column=3, pady=(0, 0))
            amino_acid2_entry = tk.Entry(table_frame, font=("Arial", 8), relief=tk.FLAT, width=15)
            amino_acid2_entry.grid(row=row + 1, column=4, pady=(0, 0))
            fraction2_entry = tk.Entry(table_frame, font=("Arial", 8), relief=tk.FLAT, width=15)
            fraction2_entry.grid(row=row + 1, column=5, pady=(0, 0))

            table_frame.rowconfigure(row + 1, pad=0)
        # create bottom frame
        bottom_frame = tk.Frame(customize_organism)
        bottom_frame.pack()

        # start_codon_label = tk.Label(bottom_frame, text="Start Codons:")
        # start_codon_label.pack(side=tk.LEFT)
        # start_codon_entry = tk.Entry(bottom_frame)
        # start_codon_entry.pack(side=tk.LEFT)
        # stop_codon_label = tk.Label(bottom_frame, text="Stop Codons:")
        # stop_codon_label.pack(side=tk.LEFT)
        # stop_codon_entry = tk.Entry(bottom_frame)
        # stop_codon_entry.pack(side=tk.LEFT)
        # check_button = tk.Button(bottom_frame, text="Check")
        # check_button.pack(side=tk.LEFT)
        save_button = tk.Button(bottom_frame, text="Save organism codon table...")
        save_button.pack(side=tk.LEFT)

        # check input
        def validate_codon_input(entry, row, col, is_amino_acid=True):
            value = entry.get()
            if value:
                if col == 1 or col == 4:  # 频率列
                    if not value.isalpha() and value != '*':
                        # self.withdraw()
                        tk.messagebox.showerror("Error", f"Invalid input in row {row + 1}, column {col + 1}. amino_acid should be a letter or *.")
                        # self.mainloop()
                        return False
                elif col == 2 or col == 5:  # 比例列
                    #if not re.match(r'^0\.\d+$', value):
                    if not re.match(r'^(0\.\d{2}|1\.00)$', value):
                        tk.messagebox.showerror("Error", f"Invalid input in row {row + 1}, column {col + 1}. fraction should be a decimal in the form 0.xx or 1.00")
                        return False
            return True

        # check table
        def check_table_content():
            errors = []
            codon_dict = {}
            for row in range(len(codons)//2):
                codon1 = codons[8 * (row // 4) + row % 4]
                codon2 = codons[8 * (row // 4) + row % 4 + 4]
                amino_acid1_entry = table_frame.grid_slaves(row=row + 1, column=1)[0]
                fraction1_entry = table_frame.grid_slaves(row=row + 1, column=2)[0]
                amino_acid2_entry = table_frame.grid_slaves(row=row + 1, column=4)[0]
                fraction2_entry = table_frame.grid_slaves(row=row + 1, column=5)[0]

                if not (validate_codon_input(amino_acid1_entry, row, 1) and validate_codon_input(fraction1_entry, row,
                                                                                                 2) and
                        validate_codon_input(amino_acid2_entry, row, 4) and validate_codon_input(fraction2_entry, row,
                                                                                                 5)):
                    continue

                amino_acid1 = amino_acid1_entry.get()
                fraction1 = fraction1_entry.get()
                amino_acid2 = amino_acid2_entry.get()
                fraction2 = fraction2_entry.get()

                if not amino_acid1 or not fraction1 or not amino_acid2 or not fraction2:
                    errors.append(f"Missing values for codons {codon1} and {codon2} in row {row + 1}.")

                # generate dictionary
                if amino_acid1 not in codon_dict:
                    codon_dict[amino_acid1] = [(codon1, fraction1)]
                else:
                    codon_dict[amino_acid1].append((codon1, fraction1))
                if amino_acid2 not in codon_dict:
                    codon_dict[amino_acid2] = [(codon2, fraction2)]
                else:
                    codon_dict[amino_acid2].append((codon2, fraction2))

            for amino_acid, codon_list in codon_dict.items():
                if amino_acid:
                    # print(f"{amino_acid}: {codon_list}")

                    total_fraction = sum([float(fraction) for _, fraction in codon_list if fraction is not None and fraction.strip()!= ""])
                    if total_fraction <= 1.01 and total_fraction >= 0.99:
                        print(f"{amino_acid}: {codon_list} (sum of fraction = {total_fraction})")
                    else:
                        print(f"{amino_acid}: {codon_list} (sum of fraction = {total_fraction})")
                        errors.append(f"{amino_acid}: {codon_list} (sum of fraction = {total_fraction})")


            # start_codon = start_codon_entry.get()
            # stop_codon = stop_codon_entry.get()
            # if not start_codon:
            #     errors.append("Missing start codon.")
            # if not stop_codon:
            #     errors.append("Missing stop codon.")

            return errors

        # save table to file
        def save_to_file():
            errors = check_table_content()
            if not errors:
                file_path = filedialog.asksaveasfilename()
                if file_path:
                    if file_path[-4:] == '.txt':
                        file_path = file_path[0:-4]
                    file = open(file_path + '.txt', 'w')
                    # file.write("organism codon table saved by Visual Codon:" + '\n')
                    for row in range(len(codons)//2):
                        codon1 = codons[8*(row // 4)+row % 4]
                        codon2 = codons[8*(row // 4)+row % 4 + 4]
                        amino_acid1_entry = table_frame.grid_slaves(row=row + 1, column=1)[0]
                        fraction1_entry = table_frame.grid_slaves(row=row + 1, column=2)[0]
                        amino_acid2_entry = table_frame.grid_slaves(row=row + 1, column=4)[0]
                        fraction2_entry = table_frame.grid_slaves(row=row + 1, column=5)[0]

                        amino_acid1 = amino_acid1_entry.get()
                        fraction1 = fraction1_entry.get()
                        amino_acid2 = amino_acid2_entry.get()
                        fraction2 = fraction2_entry.get()
                        line = codon1 + "\t" + amino_acid1 + "\t" + fraction1 + "\t"+ codon2 + "\t" + amino_acid2 + "\t"+ fraction2
                        file.write(line + '\n')

                    # data["start_codon"] = start_codon_entry.get()
                    # data["stop_codon"] = stop_codon_entry.get()

                    # file.write(self.gene_sequence + '\n')
                    file.close()

            else:
                error_message = "\n".join(errors)
                tk.messagebox.showerror("Error", error_message)

        # bind event
        for row in range(len(codons)//2):
            for col in [1, 2, 4, 5]:
                entry = table_frame.grid_slaves(row=row + 1, column=col)[0]
                entry.bind("<FocusOut>", lambda event, r=row, c=col: validate_codon_input(entry, r, c))

        # save button
        # check_button.config(command=check_table_content)
        save_button.config(command=save_to_file)

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

    Ec = """
TTT	F	0.58	TCT	S	0.17
TTC	F	0.42	TCC	S	0.15
TTA	L	0.14	TCA	S	0.14
TTG	L	0.13	TCG	S	0.14
TAT	Y	0.59	TGT	C	0.46
TAC	Y	0.41	TGC	C	0.54
TAA	*	0.61	TGA	*	0.30
TAG	*	0.09	TGG	W	1.00
CTT	L	0.12	CCT	P	0.18
CTC	L	0.10	CCC	P	0.13
CTA	L	0.04	CCA	P	0.20
CTG	L	0.47	CCG	P	0.49
CAT	H	0.57	CGT	R	0.36
CAC	H	0.43	CGC	R	0.36
CAA	Q	0.34	CGA	R	0.07
CAG	Q	0.66	CGG	R	0.11
ATT	I	0.49	ACT	T	0.19
ATC	I	0.39	ACC	T	0.4
ATA	I	0.11	ACA	T	0.17
ATG	M	1.00	ACG	T	0.25
AAT	N	0.49	AGT	S	0.16
AAC	N	0.51	AGC	S	0.25
AAA	K	0.74	AGA	R	0.07
AAG	K	0.26	AGG	R	0.04
GTT	V	0.28	GCT	A	0.18
GTC	V	0.20	GCC	A	0.26
GTA	V	0.17	GCA	A	0.23
GTG	V	0.35	GCG	A	0.33
GAT	D	0.63	GGT	G	0.35
GAC	D	0.37	GGC	G	0.37
GAA	E	0.68	GGA	G	0.13
GAG	E	0.32	GGG	G	0.15
"""

    Yeast = """
TTT	F	0.59	TCT	S	0.26
TTC	F	0.41	TCC	S	0.16
TTA	L	0.28	TCA	S	0.21
TTG	L	0.29	TCG	S	0.10
TAT	Y	0.56	TGT	C	0.63
TAC	Y	0.44	TGC	C	0.37
TAA	*	0.48	TGA	*	0.29
TAG	*	0.24	TGG	W	1.00
CTT	L	0.13	CCT	P	0.31
CTC	L	0.06	CCC	P	0.15
CTA	L	0.14	CCA	P	0.41
CTG	L	0.11	CCG	P	0.12
CAT	H	0.64	CGT	R	0.15
CAC	H	0.36	CGC	R	0.06
CAA	Q	0.69	CGA	R	0.07
CAG	Q	0.31	CGG	R	0.04
ATT	I	0.46	ACT	T	0.35
ATC	I	0.26	ACC	T	0.22
ATA	I	0.27	ACA	T	0.30
ATG	M	1.00	ACG	T	0.13
AAT	N	0.59	AGT	S	0.16
AAC	N	0.41	AGC	S	0.11
AAA	K	0.58	AGA	R	0.48
AAG	K	0.42	AGG	R	0.21
GTT	V	0.39	GCT	A	0.38
GTC	V	0.21	GCC	A	0.22
GTA	V	0.21	GCA	A	0.29
GTG	V	0.19	GCG	A	0.11
GAT	D	0.65	GGT	G	0.47
GAC	D	0.35	GGC	G	0.19
GAA	E	0.71	GGA	G	0.22
GAG	E	0.29	GGG	G	0.12
"""

    Insect = """
TTT	F	0.25	TCT	S	0.19
TTC	F	0.75	TCC	S	0.24
TTA	L	0.07	TCA	S	0.16
TTG	L	0.20	TCG	S	0.12
TAT	Y	0.25	TGT	C	0.35
TAC	Y	0.75	TGC	C	0.65
TAA	*	0.68	TGA	*	0.17
TAG	*	0.15	TGG	W	1.00
CTT	L	0.13	CCT	P	0.29
CTC	L	0.22	CCC	P	0.31
CTA	L	0.08	CCA	P	0.24
CTG	L	0.31	CCG	P	0.16
CAT	H	0.32	CGT	R	0.25
CAC	H	0.68	CGC	R	0.24
CAA	Q	0.39	CGA	R	0.08
CAG	Q	0.61	CGG	R	0.05
ATT	I	0.29	ACT	T	0.27
ATC	I	0.59	ACC	T	0.36
ATA	I	0.12	ACA	T	0.21
ATG	M	1.00	ACG	T	0.16
AAT	N	0.29	AGT	S	0.12
AAC	N	0.71	AGC	S	0.17
AAA	K	0.31	AGA	R	0.16
AAG	K	0.69	AGG	R	0.21
GTT	V	0.20	GCT	A	0.36
GTC	V	0.30	GCC	A	0.30
GTA	V	0.15	GCA	A	0.16
GTG	V	0.35	GCG	A	0.18
GAT	D	0.37	GGT	G	0.36
GAC	D	0.63	GGC	G	0.32
GAA	E	0.41	GGA	G	0.27
GAG	E	0.59	GGG	G	0.05
"""

    Cele="""
TTT	F	0.50	TCT	S	0.21
TTC	F	0.50	TCC	S	0.13
TTA	L	0.12	TCA	S	0.25
TTG	L	0.23	TCG	S	0.15
TAT	Y	0.56	TGT	C	0.55
TAC	Y	0.44	TGC	C	0.45
TAA	*	0.44	TGA	*	0.39
TAG	*	0.17	TGG	W	1.00
CTT	L	0.24	CCT	P	0.18
CTC	L	0.17	CCC	P	0.09
CTA	L	0.09	CCA	P	0.53
CTG	L	0.14	CCG	P	0.20
CAT	H	0.61	CGT	R	0.21
CAC	H	0.39	CGC	R	0.10
CAA	Q	0.66	CGA	R	0.23
CAG	Q	0.34	CGG	R	0.09
ATT	I	0.53	ACT	T	0.33
ATC	I	0.31	ACC	T	0.18
ATA	I	0.16	ACA	T	0.34
ATG	M	1.00	ACG	T	0.15
AAT	N	0.62	AGT	S	0.15
AAC	N	0.38	AGC	S	0.10
AAA	K	0.59	AGA	R	0.29
AAG	K	0.41	AGG	R	0.08
GTT	V	0.39	GCT	A	0.36
GTC	V	0.22	GCC	A	0.20
GTA	V	0.16	GCA	A	0.31
GTG	V	0.23	GCG	A	0.13
GAT	D	0.68	GGT	G	0.20
GAC	D	0.33	GGC	G	0.12
GAA	E	0.62	GGA	G	0.59
GAG	E	0.38	GGG	G	0.08
"""

    Droso="""
TTT	F	0.37	TCT	S	0.08
TTC	F	0.63	TCC	S	0.24
TTA	L	0.05	TCA	S	0.09
TTG	L	0.18	TCG	S	0.20
TAT	Y	0.37	TGT	C	0.29
TAC	Y	0.63	TGC	C	0.71
TAA	*	0.42	TGA	*	0.26
TAG	*	0.32	TGG	W	1.00
CTT	L	0.10	CCT	P	0.13
CTC	L	0.15	CCC	P	0.33
CTA	L	0.09	CCA	P	0.25
CTG	L	0.43	CCG	P	0.29
CAT	H	0.40	CGT	R	0.16
CAC	H	0.60	CGC	R	0.33
CAA	Q	0.30	CGA	R	0.15
CAG	Q	0.70	CGG	R	0.15
ATT	I	0.34	ACT	T	0.17
ATC	I	0.47	ACC	T	0.38
ATA	I	0.19	ACA	T	0.19
ATG	M	1.00	ACG	T	0.26
AAT	N	0.44	AGT	S	0.14
AAC	N	0.56	AGC	S	0.25
AAA	K	0.29	AGA	R	0.09
AAG	K	0.71	AGG	R	0.11
GTT	V	0.18	GCT	A	0.19
GTC	V	0.24	GCC	A	0.45
GTA	V	0.11	GCA	A	0.17
GTG	V	0.47	GCG	A	0.19
GAT	D	0.53	GGT	G	0.21
GAC	D	0.47	GGC	G	0.43
GAA	E	0.33	GGA	G	0.29
GAG	E	0.67	GGG	G	0.07
"""

    Human = """
TTT	F	0.45	TCT	S	0.18
TTC	F	0.55	TCC	S	0.22
TTA	L	0.07	TCA	S	0.15
TTG	L	0.13	TCG	S	0.06
TAT	Y	0.43	TGT	C	0.45
TAC	Y	0.57	TGC	C	0.55
TAA	*	0.28	TGA	*	0.52
TAG	*	0.20	TGG	W	1.00
CTT	L	0.13	CCT	P	0.28
CTC	L	0.20	CCC	P	0.33
CTA	L	0.07	CCA	P	0.27
CTG	L	0.41	CCG	P	0.11
CAT	H	0.41	CGT	R	0.08
CAC	H	0.59	CGC	R	0.19
CAA	Q	0.25	CGA	R	0.11
CAG	Q	0.75	CGG	R	0.21
ATT	I	0.36	ACT	T	0.24
ATC	I	0.48	ACC	T	0.36
ATA	I	0.16	ACA	T	0.28
ATG	M	1.00	ACG	T	0.12
AAT	N	0.46	AGT	S	0.15
AAC	N	0.54	AGC	S	0.24
AAA	K	0.42	AGA	R	0.20
AAG	K	0.58	AGG	R	0.20
GTT	V	0.18	GCT	A	0.26
GTC	V	0.24	GCC	A	0.40
GTA	V	0.11	GCA	A	0.23
GTG	V	0.47	GCG	A	0.11
GAT	D	0.46	GGT	G	0.16
GAC	D	0.54	GGC	G	0.34
GAA	E	0.42	GGA	G	0.25
GAG	E	0.58	GGG	G	0.25
"""

    Mouse = """
TTT	F	0.43	TCT	S	0.19
TTC	F	0.57	TCC	S	0.22
TTA	L	0.06	TCA	S	0.14
TTG	L	0.13	TCG	S	0.05
TAT	Y	0.43	TGT	C	0.48
TAC	Y	0.58	TGC	C	0.52
TAA	*	0.26	TGA	*	0.52
TAG	*	0.22	TGG	W	1.00
CTT	L	0.13	CCT	P	0.30
CTC	L	0.20	CCC	P	0.31
CTA	L	0.08	CCA	P	0.28
CTG	L	0.39	CCG	P	0.10
CAT	H	0.40	CGT	R	0.09
CAC	H	0.60	CGC	R	0.18
CAA	Q	0.25	CGA	R	0.12
CAG	Q	0.75	CGG	R	0.19
ATT	I	0.34	ACT	T	0.25
ATC	I	0.50	ACC	T	0.35
ATA	I	0.16	ACA	T	0.29
ATG	M	1.00	ACG	T	0.11
AAT	N	0.43	AGT	S	0.15
AAC	N	0.57	AGC	S	0.24
AAA	K	0.39	AGA	R	0.21
AAG	K	0.61	AGG	R	0.22
GTT	V	0.17	GCT	A	0.29
GTC	V	0.25	GCC	A	0.38
GTA	V	0.12	GCA	A	0.23
GTG	V	0.46	GCG	A	0.10
GAT	D	0.44	GGT	G	0.18
GAC	D	0.56	GGC	G	0.33
GAA	E	0.40	GGA	G	0.26
GAG	E	0.60	GGG	G	0.23
"""

    Rat = """
TTT	F	0.41	TCT	S	0.18
TTC	F	0.59	TCC	S	0.23
TTA	L	0.06	TCA	S	0.14
TTG	L	0.12	TCG	S	0.06
TAT	Y	0.40	TGT	C	0.45
TAC	Y	0.60	TGC	C	0.56
TAA	*	0.27	TGA	*	0.50
TAG	*	0.23	TGG	W	1.00
CTT	L	0.12	CCT	P	0.30
CTC	L	0.21	CCC	P	0.32
CTA	L	0.07	CCA	P	0.27
CTG	L	0.42	CCG	P	0.11
CAT	H	0.38	CGT	R	0.09
CAC	H	0.62	CGC	R	0.18
CAA	Q	0.25	CGA	R	0.12
CAG	Q	0.75	CGG	R	0.20
ATT	I	0.32	ACT	T	0.23
ATC	I	0.54	ACC	T	0.37
ATA	I	0.14	ACA	T	0.28
ATG	M	1.00	ACG	T	0.12
AAT	N	0.40	AGT	S	0.15
AAC	N	0.60	AGC	S	0.25
AAA	K	0.37	AGA	R	0.20
AAG	K	0.63	AGG	R	0.21
GTT	V	0.16	GCT	A	0.28
GTC	V	0.26	GCC	A	0.40
GTA	V	0.11	GCA	A	0.22
GTG	V	0.47	GCG	A	0.10
GAT	D	0.42	GGT	G	0.17
GAC	D	0.58	GGC	G	0.34
GAA	E	0.39	GGA	G	0.25
GAG	E	0.61	GGG	G	0.24
"""

    Pig = """
TTT	F	0.38	TCT	S	0.15
TTC	F	0.62	TCC	S	0.26
TTA	L	0.06	TCA	S	0.15
TTG	L	0.10	TCG	S	0.06
TAT	Y	0.35	TGT	C	0.39
TAC	Y	0.65	TGC	C	0.61
TAA	*	0.13	TGA	*	0.79
TAG	*	0.08	TGG	W	1.00
CTT	L	0.10	CCT	P	0.24
CTC	L	0.21	CCC	P	0.35
CTA	L	0.13	CCA	P	0.27
CTG	L	0.40	CCG	P	0.13
CAT	H	0.34	CGT	R	0.07
CAC	H	0.66	CGC	R	0.22
CAA	Q	0.25	CGA	R	0.12
CAG	Q	0.75	CGG	R	0.20
ATT	I	0.30	ACT	T	0.19
ATC	I	0.53	ACC	T	0.41
ATA	I	0.18	ACA	T	0.26
ATG	M	1.00	ACG	T	0.13
AAT	N	0.37	AGT	S	0.12
AAC	N	0.63	AGC	S	0.26
AAA	K	0.40	AGA	R	0.19
AAG	K	0.60	AGG	R	0.20
GTT	V	0.14	GCT	A	0.24
GTC	V	0.27	GCC	A	0.45
GTA	V	0.12	GCA	A	0.20
GTG	V	0.48	GCG	A	0.11
GAT	D	0.38	GGT	G	0.14
GAC	D	0.62	GGC	G	0.36
GAA	E	0.38	GGA	G	0.26
GAG	E	0.62	GGG	G	0.24
"""

    Pichia = """
TTT	F	0.56	TCT	S	0.29
TTC	F	0.44	TCC	S	0.20
TTA	L	0.15	TCA	S	0.19
TTG	L	0.33	TCG	S	0.09
TAT	Y	0.45	TGT	C	0.65
TAC	Y	0.55	TGC	C	0.35
TAA	*	0.53	TGA	*	0.18
TAG	*	0.29	TGG	W	1.00
CTT	L	0.17	CCT	P	0.35
CTC	L	0.08	CCC	P	0.16
CTA	L	0.12	CCA	P	0.40
CTG	L	0.16	CCG	P	0.10
CAT	H	0.54	CGT	R	0.16
CAC	H	0.46	CGC	R	0.05
CAA	Q	0.62	CGA	R	0.11
CAG	Q	0.38	CGG	R	0.05
ATT	I	0.51	ACT	T	0.40
ATC	I	0.31	ACC	T	0.24
ATA	I	0.18	ACA	T	0.25
ATG	M	1.00	ACG	T	0.11
AAT	N	0.48	AGT	S	0.15
AAC	N	0.52	AGC	S	0.09
AAA	K	0.47	AGA	R	0.47
AAG	K	0.53	AGG	R	0.16
GTT	V	0.42	GCT	A	0.45
GTC	V	0.23	GCC	A	0.25
GTA	V	0.16	GCA	A	0.24
GTG	V	0.20	GCG	A	0.06
GAT	D	0.59	GGT	G	0.43
GAC	D	0.41	GGC	G	0.14
GAA	E	0.58	GGA	G	0.32
GAG	E	0.42	GGG	G	0.10
"""

    Arabid = """
TTT	F	0.51	TCT	S	0.28
TTC	F	0.49	TCC	S	0.13
TTA	L	0.13	TCA	S	0.20
TTG	L	0.22	TCG	S	0.10
TAT	Y	0.52	TGT	C	0.59
TAC	Y	0.48	TGC	C	0.41
TAA	*	0.36	TGA	*	0.44
TAG	*	0.20	TGG	W	1.00
CTT	L	0.26	CCT	P	0.38
CTC	L	0.17	CCC	P	0.11
CTA	L	0.11	CCA	P	0.33
CTG	L	0.11	CCG	P	0.17
CAT	H	0.61	CGT	R	0.17
CAC	H	0.39	CGC	R	0.07
CAA	Q	0.56	CGA	R	0.12
CAG	Q	0.44	CGG	R	0.09
ATT	I	0.41	ACT	T	0.34
ATC	I	0.35	ACC	T	0.20
ATA	I	0.24	ACA	T	0.30
ATG	M	1.00	ACG	T	0.15
AAT	N	0.52	AGT	S	0.16
AAC	N	0.48	AGC	S	0.13
AAA	K	0.48	AGA	R	0.35
AAG	K	0.52	AGG	R	0.20
GTT	V	0.41	GCT	A	0.44
GTC	V	0.19	GCC	A	0.16
GTA	V	0.15	GCA	A	0.27
GTG	V	0.26	GCG	A	0.14
GAT	D	0.68	GGT	G	0.34
GAC	D	0.32	GGC	G	0.14
GAA	E	0.52	GGA	G	0.37
GAG	E	0.49	GGG	G	0.15
"""

    Strep = """
TTT	F	0.02	TCT	S	0.01
TTC	F	0.98	TCC	S	0.41
TTA	L	0.00	TCA	S	0.02
TTG	L	0.02	TCG	S	0.28
TAT	Y	0.05	TGT	C	0.09
TAC	Y	0.95	TGC	C	0.91
TAA	*	0.03	TGA	*	0.80
TAG	*	0.17	TGG	W	1.00
CTT	L	0.02	CCT	P	0.02
CTC	L	0.36	CCC	P	0.41
CTA	L	0.00	CCA	P	0.02
CTG	L	0.60	CCG	P	0.54
CAT	H	0.07	CGT	R	0.07
CAC	H	0.93	CGC	R	0.47
CAA	Q	0.05	CGA	R	0.03
CAG	Q	0.95	CGG	R	0.38
ATT	I	0.02	ACT	T	0.02
ATC	I	0.96	ACC	T	0.65
ATA	I	0.02	ACA	T	0.03
ATG	M	1.00	ACG	T	0.31
AAT	N	0.04	AGT	S	0.03
AAC	N	0.96	AGC	S	0.25
AAA	K	0.05	AGA	R	0.01
AAG	K	0.95	AGG	R	0.04
GTT	V	0.02	GCT	A	0.02
GTC	V	0.55	GCC	A	0.57
GTA	V	0.03	GCA	A	0.04
GTG	V	0.41	GCG	A	0.36
GAT	D	0.05	GGT	G	0.10
GAC	D	0.95	GGC	G	0.64
GAA	E	0.15	GGA	G	0.08
GAG	E	0.85	GGG	G	0.19
"""

    Zea = """
TTT	F	0.37	TCT	S	0.17
TTC	F	0.63	TCC	S	0.22
TTA	L	0.08	TCA	S	0.15
TTG	L	0.15	TCG	S	0.14
TAT	Y	0.37	TGT	C	0.34
TAC	Y	0.63	TGC	C	0.66
TAA	*	0.24	TGA	*	0.44
TAG	*	0.32	TGG	W	1.00
CTT	L	0.18	CCT	P	0.24
CTC	L	0.25	CCC	P	0.24
CTA	L	0.08	CCA	P	0.26
CTG	L	0.25	CCG	P	0.26
CAT	H	0.43	CGT	R	0.12
CAC	H	0.57	CGC	R	0.23
CAA	Q	0.39	CGA	R	0.09
CAG	Q	0.61	CGG	R	0.15
ATT	I	0.33	ACT	T	0.24
ATC	I	0.47	ACC	T	0.33
ATA	I	0.20	ACA	T	0.22
ATG	M	1.00	ACG	T	0.21
AAT	N	0.40	AGT	S	0.11
AAC	N	0.60	AGC	S	0.21
AAA	K	0.30	AGA	R	0.16
AAG	K	0.70	AGG	R	0.25
GTT	V	0.24	GCT	A	0.25
GTC	V	0.29	GCC	A	0.33
GTA	V	0.11	GCA	A	0.19
GTG	V	0.36	GCG	A	0.23
GAT	D	0.44	GGT	G	0.21
GAC	D	0.56	GGC	G	0.39
GAA	E	0.36	GGA	G	0.20
GAG	E	0.64	GGG	G	0.21
"""

    Nicotia = """
TTT	F	0.58	TCT	S	0.26
TTC	F	0.42	TCC	S	0.14
TTA	L	0.14	TCA	S	0.23
TTG	L	0.24	TCG	S	0.07
TAT	Y	0.57	TGT	C	0.57
TAC	Y	0.43	TGC	C	0.43
TAA	*	0.41	TGA	*	0.41
TAG	*	0.19	TGG	W	1.00
CTT	L	0.26	CCT	P	0.37
CTC	L	0.14	CCC	P	0.13
CTA	L	0.10	CCA	P	0.40
CTG	L	0.12	CCG	P	0.10
CAT	H	0.61	CGT	R	0.15
CAC	H	0.39	CGC	R	0.08
CAA	Q	0.58	CGA	R	0.11
CAG	Q	0.42	CGG	R	0.08
ATT	I	0.50	ACT	T	0.39
ATC	I	0.25	ACC	T	0.19
ATA	I	0.25	ACA	T	0.33
ATG	M	1.00	ACG	T	0.09
AAT	N	0.60	AGT	S	0.17
AAC	N	0.40	AGC	S	0.13
AAA	K	0.49	AGA	R	0.32
AAG	K	0.51	AGG	R	0.26
GTT	V	0.41	GCT	A	0.44
GTC	V	0.17	GCC	A	0.17
GTA	V	0.17	GCA	A	0.31
GTG	V	0.25	GCG	A	0.08
GAT	D	0.68	GGT	G	0.34
GAC	D	0.32	GGC	G	0.17
GAA	E	0.55	GGA	G	0.34
GAG	E	0.45	GGG	G	0.15
"""

    Saccha = """
TTT	F	0.59	TCT	S	0.26
TTC	F	0.41	TCC	S	0.16
TTA	L	0.28	TCA	S	0.21
TTG	L	0.29	TCG	S	0.10
TAT	Y	0.56	TGT	C	0.63
TAC	Y	0.44	TGC	C	0.37
TAA	*	0.48	TGA	*	0.30
TAG	*	0.22	TGG	W	1.00
CTT	L	0.13	CCT	P	0.31
CTC	L	0.06	CCC	P	0.15
CTA	L	0.14	CCA	P	0.42
CTG	L	0.11	CCG	P	0.12
CAT	H	0.64	CGT	R	0.14
CAC	H	0.36	CGC	R	0.06
CAA	Q	0.69	CGA	R	0.07
CAG	Q	0.31	CGG	R	0.04
ATT	I	0.46	ACT	T	0.35
ATC	I	0.26	ACC	T	0.22
ATA	I	0.27	ACA	T	0.30
ATG	M	1.00	ACG	T	0.14
AAT	N	0.59	AGT	S	0.16
AAC	N	0.41	AGC	S	0.11
AAA	K	0.58	AGA	R	0.48
AAG	K	0.42	AGG	R	0.21
GTT	V	0.39	GCT	A	0.38
GTC	V	0.21	GCC	A	0.22
GTA	V	0.21	GCA	A	0.29
GTG	V	0.19	GCG	A	0.11
GAT	D	0.65	GGT	G	0.47
GAC	D	0.35	GGC	G	0.19
GAA	E	0.70	GGA	G	0.22
GAG	E	0.30	GGG	G	0.12
"""

    Cricet = """
TTT	F	0.47	TCT	S	0.22
TTC	F	0.53	TCC	S	0.22
TTA	L	0.07	TCA	S	0.14
TTG	L	0.14	TCG	S	0.05
TAT	Y	0.44	TGT	C	0.47
TAC	Y	0.56	TGC	C	0.53
TAA	*	0.26	TGA	*	0.53
TAG	*	0.22	TGG	W	1.00
CTT	L	0.13	CCT	P	0.31
CTC	L	0.19	CCC	P	0.32
CTA	L	0.08	CCA	P	0.29
CTG	L	0.39	CCG	P	0.08
CAT	H	0.44	CGT	R	0.11
CAC	H	0.56	CGC	R	0.18
CAA	Q	0.24	CGA	R	0.14
CAG	Q	0.76	CGG	R	0.19
ATT	I	0.35	ACT	T	0.26
ATC	I	0.51	ACC	T	0.37
ATA	I	0.14	ACA	T	0.29
ATG	M	1.00	ACG	T	0.08
AAT	N	0.45	AGT	S	0.15
AAC	N	0.55	AGC	S	0.22
AAA	K	0.39	AGA	R	0.19
AAG	K	0.61	AGG	R	0.19
GTT	V	0.18	GCT	A	0.32
GTC	V	0.24	GCC	A	0.37
GTA	V	0.12	GCA	A	0.23
GTG	V	0.46	GCG	A	0.07
GAT	D	0.47	GGT	G	0.20
GAC	D	0.53	GGC	G	0.34
GAA	E	0.41	GGA	G	0.25
GAG	E	0.59	GGG	G	0.21
"""

if __name__ == "__main__":
    # root = tk.Tk()
    app = App()
    app.mainloop()

