

### Visual codon: A user-friendly Python program for viewing and optimizing gene GC content

This program is used to display, modify, and compare the codons of gene. It shows the codons and GC content of the gene in the form of visual graphics, and the user can modify the codons to optimize GC content.

At the same time, the program also supports the codon comparison function of two gene sequences coding for the same protein. When using the comparison function, the program will enter the read-only mode, and the codon cannot be modified. If you need to modify it, please reopen the sequence file through the menu 'Open Gene 1 to edit'.

The program has built-in support for 16 species including Escherichia coil, Yeast, Insect, C. elegans, Drosophila melanogaster, Human, Mouse, Rat, Pig, Pichia pastoris, Arabidopsis thaliana, Streptomyces, Zea mays (Maize), Nicotiana tabacum (Tabacco), Saccharomyces cerevisiae (gbpln), and Cricetulus griseus (CHO). Please visit the GenScript website for detail at https://www.genscript.com/tools/codon-frequency-table.

The program also supports custom host species. Users can configure the codon table of the host species  through the program interface, or manually edit the txt file to achieve it.
This program can be run on macOS/Linux/Windows and the results are the same. Please follow the instructions below to install and run the program.

The check sequence function of this program calls the EMBOSS software (https://emboss.sourceforge.net/). If you need to use the check sequence function, please make sure that the EMBOSS software is installed in advance. The EMBOSS version used by the author is 6.6.0.0. Please refer to https://emboss.sourceforge.net/apps/ for installation instructions.




### I. Installation
#### <font color=red> in macOS/Linux </font>
(1) make sure you have installed Python 3.12.
Note: The author's environment is Python 3.12, other versions of Python and matplotlib may also be able to run this program.

```zsh
% python3.12
Python 3.12.1 (v3.12.1:2305ca5144, Dec  7 2023, 17:23:38) [Clang 13.0.0 (clang-1300.0.29.30)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
```

Note that some computers may have more than one Python version installed, you need to check the version number of the python command. If you enter 'python3.12' command prompt and do not find it, you can try to use the 'python' command directly.

The following steps are based on 'python 3.12' command, and you can decide whether to change it to 'python' command based on your specific environment.

(2) make sure you have installed the dependent package matplotlib 3.9.2.

```zsh
% python3.12 -m pip list | grep -E '^(matplot)'
matplotlib                                3.9.2
```

The above output indicates that the matplotlib 3.9.2 package has been installed.
If the dependency is not installed, install it with the following command:

```zsh
% python3.12 -m pip install matplotlib==3.9.2
```

If you have installed the matplotlib package but the version is inconsistent, please use the following command to update the version of matplotlib:

```zsh
% python3.12 -m pip install --upgrade matplotlib==3.9.2
```

Confirm whether the version is correct after the update is complete.

```zsh
% python3.12 -m pip show matplotlib
Name: matplotlib
Version: 3.9.2
Summary: Python plotting package
Home-page:
Author: John D. Hunter, Michael Droettboom
Author-email: Unknown <matplotlib-users@python.org>
License: License agreement for matplotlib versions 1.3.0 and later
```

In general, pip will automatically install matplotlib dependencies (Requires: contourpy, cycler, fonttools, kiwisolver, numpy, packaging, pillow, pyparsing, python-dateutil) when installing matplotlib.

If pip does not automatically install them in your environment, you need to install these dependencies manually.

(3) Use git clone or Download ZIP to download program file 'visual_codon.py' and test sequence files to your local computer.

GitHub address: https://github.com/shiqiang-lin/visual_codon

(4) Save the 'visual_codon.py' to a directory with read and write permissions. For example:

```zsh
% cd ~            # enter Home directory
% cd Desktop      # enter Desktop
% mkdir codons    # create a folder 'codons'. Other names can be used too.
```

Copy the downloaded 'visual_codon.py' to the 'codons' directory. This can be done by dragging or right-clicking with the mouse.

#### <font color=red> in Windows </font>
(1) make sure you have installed Python 3.12. 
Note: The author's environment is Python 3.12, other versions of Python and matplotlib may also be able to run this program.

```zsh
> python3.12
Python 3.12.1 (v3.12.1:2305ca5144, Dec  7 2023, 17:23:38) [MSC v.1940 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
```

Note that some computers may have more than one Python version installed, you need to check the version number of the python command. If you enter 'python3.12' command prompt and do not find it, you can try to use the 'python' command directly.

The following steps are based on 'python 3.12' command, and you can decide whether to change it to 'python' command based on your specific environment.

(2) make sure you have installed the dependent package matplotlib 3.9.2.

```zsh
> python3.12 -m pip list | findstr "matplot"
matplotlib 3.9.2      
```

The above output indicates that the matplotlib 3.9.2 package has been installed. If the dependency is not installed, install it with the following command:

```zsh
> python3.12 -m pip install matplotlib==3.9.2
```

If you have installed the matplotlib package but the version is inconsistent, please use the following command to update the version of matplotlib:

```zsh
> python3.12 -m pip install --upgrade matplotlib==3.9.2
```

Confirm whether the version is correct after the update is complete.

```zsh
% python3.12 -m pip show matplotlib
Name: matplotlib
Version: 3.9.2
Summary: Python plotting package
Home-page:
Author: John D. Hunter, Michael Droettboom
Author-email: Unknown <matplotlib-users@python.org>
License: License agreement for matplotlib versions 1.3.0 and later
```
In general, pip will automatically install matplotlib dependencies (Requires: contourpy, cycler, fonttools, kiwisolver, numpy, packaging, pillow, pyparsing, python-dateutil) when installing matplotlib.

If pip does not automatically install them in your environment, you need to install these dependencies manually.

(3) Use git clone or Download ZIP to download program file 'visual_codon.py' and test sequence files to your local computer.

GitHub address: https://github.com/shiqiang-lin/visual_codon

(4) Save the 'visual_codon.py' to a directory with read and write permissions. For example:
Right-click on the desktop to create a new folder 'codons', copy the downloaded 'visual_codon.py' to 'codons', and you can use the mouse to drag or right-click to copy and paste.

### II. Usage

Here are the steps to start the program:
#### <font color=red> in macOS/Linux </font>

Operate the terminal, go to the Desktop/codons directory, and start the program.

```zsh
% cd ~
% cd Desktop/codons
% python3.12 visual_codon.py
```

#### <font color=red> in Windows </font>
Open the Command Prompt window, cmd.exe, go to the C:/Users/44139/Desktop/codons directory (this directory is the author's environment, please modify it according to the actual path), and start the program.

```zsh
> cd C:/Users/44139/Desktop/codons
> python3.12 visual_codon.py
```

#### <font color=red> in X11 forwarding </font>

This program supports X11 forwarding to connect to remote Linux servers. Here we show how to use a Windows 11 client to connect to a remote Linux server using the SecureCRT client.


(1) Install the Xming program on Windows 11.

(2) Use SecureCRT to log in to the remote Linux server with SSH2.

(3) After logging in to the Linux server, set the environment variable DISPLAY. The 192.168.31.129 is the IP address of the Windows computer of the author, please replace it with the IP address of your own computer.

```zsh
% export DISPLAY=192.168.31.129:0.0
% 
```
(4) Run the xclock program to see if the X11 forwarding configuration is normal. If the xclock clock can be displayed, the configuration is correct.

```zsh
% xclock
% 
```
(5) Go to the Desktop/codons directory, and launch the program.

```zsh
% cd ~
% cd Desktop/codons
% python3.12 visual_codon.py
```

After launching the program, the program will open the Select a host organism interface, and here are the instructions for using it.

#### <font color=red> in macOS/Linux/Windows </font>

(1) On the host species selection page, you can select the following 16 species from the drop-down box, "Escherichia coil", "Yeast", "Insect", "C. elegans", "Drosophila melanogaster", "Human", "Mouse", "Rat", "Pig", "Pichia pastoris", "Arabidopsis thaliana", "Streptomyces", "Zea mays (Maize)", "Nicotiana tabacum (Tabacco)", "Saccharomyces cerevisiae (gbpln)", "Cricetulus griseus (CHO)" (Please refer to https://www.genscript.com/tools/codon-frequency-table for detail), or click Open File... to open the saved host species file.

(2) Click File->'Open Gene 1' to edit to open the gene sequence.
The file name suffix should be .fasta or .fa. After opening the file, the program will check the validity of the sequence file, including the length must be greater than or equal to 60 and be a multiple of 3.
There should be no stop codon in the middle of the sequence, and the sequence cannot contain characters other than ATCG. The program supports lowercase characters for sequence files, which will be automatically converted to uppercase when processed by the program.

Note: When opening the gene file, the program will pop up a prompt message about the host species used, to help users confirm. The host species will also be displayed at the end of the program GUI title after opening the gene file.
If you need to change the host species, close the program and select the host species when you reopen the program.

(3) Graphical display.
After the sequence file is opened, the TreeView on the main interface of the program will display each amino acid and codon information. Before the codon is modified, the original codon, original rank, and original GC% are the same as the optimized codon, optimized rank, and optimized GC%.

The selected amino acid and codon information are displayed in the middle of the screen. RadioButtons display all the optional codons of the amino acid for the user to choose to modify.

At the bottom of the interface is the GC content plot, the abscissa is codon No., and the ordinate is GC%. When you click the left mouse button on the graphic line, the treeview will automatically jump to the corresponding position.

At the bottom left of the graph is the navigation toolbar of the graph, through which you can drag, zoom in, and save the graph. For more information on how to navigate the toolbar, please refer to the toolbar manual at
https://matplotlib.org/stable/users/explain/figure/interactive.html#interactive-navigation

(4) Modify codons.
Select the amino acid you want to modify by clicking on the GC content graph or by clicking on a row in the treeview. In the radiobutton in the middle of the interface, select the codon you want to change, and click Update codon to modify.

After the modification is completed, optimized codon, optimized rank, optimized GC% will be updated to the new values. At the same time, the graphic lines at the bottom of the interface will also be refreshed, with the brown line being Gene 1 (original gene) and the blue line being the modified gene (Gene 1 optimized).

Users can modify multiple codons to determine the effect of the update based on the shape of the graphic lines.

(5) Click File->'Save optimized gene' to save the optimized gene.

After the modification of the codon is completed, click File->'Save optimized gene' to save the optimized gene, and the new gene will be saved as a fasta file.

(6) Click File->'Export table to txt' to save the treeview list.

By clicking File->'Export table to txt', you can save the treeview list as a txt file, and the txt file can be imported to excel (Windows) or Numbers (macOS) for statistics and other operations.

(7) Click File->'Export changed codons' to save the updated rows in the treeview list.

By clicking File->'Export changed codons', you can save the updated rows in the treeview list as a txt file, and the txt file can be imported to excel (Windows) or Numbers (macOS) for statistics and other operations.

(8) Click File->'Import Gene 2 to compare' and compare the two gene sequences.

If you modify an open Gene 1, the File->'Import Gene 2' to compare menu item will be disabled.
If you want to see the difference between the GC% of the two gene sequences, you can click File->> 'Import Gene 2 to compare' after File- 'Open Gene 1 to edit' to import Gene 2.

When Gene 2 is introduced, it is also checked for validity, and Gene 1 and Gene 2 must be the same length. The graph at the bottom of the interface will show two GC content lines. The brown line represents Gene 1, and the blue line represents Gene 2. Users can also operate the navigation toolbar to drag, zoom in, zoom in, and save graphics. It is worth noting that when comparing two genes, the program will enter comparision mode [ReadOnly],
The codon cannot be modified in this mode (the Update codon button is disabled).

(9) Click Edit->'Adjust graph font' to modify the font size and interval of the X axis
The program automatically adjusts the interval according to the length of the gene sequence, if the user still wants to adjust it can be done through this menu item.

(10) Click Edit->'Customize organism' to customize the host species codon table, and the user can fill in the form with the amino acids and fractions corresponding to the codons.
When you finish filling in the form and save it to a txt file, the program will check whether the user's filling in content meets the requirements.

The judging conditions include: the amino acid must be a letter or * (stop codon), the fraction must be 0.xx or 1.00, and the sum of the fractions of each codon corresponding to the same amino acid is required to be in the range of [0.99, 1.01].

Note: The visual_codon program directory provides the 'example_of_customized_codon_table.txt', and users can also modify this sample file directly.

After the modification is completed, open the visual codon program, click Open File on the host species selection screen... Select the host species file that you just modified.

Note: When modifying the example_of_customized_codon_table.txt, you must ensure that the amino acids, codons, and fractions filled in are correct.

(11) Click Edit->'Check sequence' to verify the DNA sequence.
The Check sequence interface provides five EMBOSS commands: equicktandem, etandem, palindrome, einverted, and restrict.(https://emboss.sourceforge.net/).

PLEASE MAKE SURE THAT EMBOSS IS INSTALLED CORRECTLY BEFORE USING THIS FEATURE.

(12) Click File->'Reinitialize' to reset the gene file, after which the new gene file can be opened, and the host species has not been changed.

If you need to change the host species, please close the program and select the desired host species on the host species selection screen when you reopen the program.
