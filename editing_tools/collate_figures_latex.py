#!/usr/bin/python
# --------------------------------------------------------
# Utility to stitch a pdf document of figures and captions
# from a list of figure files and caption files,
# with the main purpose of stitching together supplements
# Carles Boix 09/17/19
# --------------------------------------------------------
import os  # For system + paths
import re  # For captions
import fire  # To make cli tool


class collate_figures_pdf(object):
    """Utility to stitching together a PDF from a list of figures and captions.

    Run as:
        python collate_figures_latex.py main \
            --figurelist FIGURELIST.tsv --outfile OUT.pdf

    Args:
        figurelist (str):
            Tab-separated table with information about figures:
            Column 1: Figure filenames
            Column 2: Figure captions, in latex style.
            Columns 3+: Can be anything, will not read these.
        outfile (str):
            Either the intended pdf or tex file. Will strip extension.

    Note:
        Run this in the figure directory. Otherwise, pdflatex
        will not find figures unless they contain a full path.
    """
    def __init__(self, figurelist, outfile, bold_captions=False):
        self.figurelist = figurelist
        # NOTE: outfile can be tex or pdf, will strip extension.
        self.outprefix = os.path.splitext(outfile)[0]
        self.outtex = self.outprefix + ".tex"
        self.outpdf = self.outprefix + ".pdf"
        # If want to auto modify captions:
        self.bold_captions = bold_captions

    def main(self):
        self.read_lists()
        if self.bold_captions:
            self.modify_captions()
        self.write_tex_file()
        self.compile_tex_file()

    def read_lists(self):
        self.figures = []
        self.captions = []
        with open(self.figurelist, 'r') as f:
            for line in f:
                line = line.strip("\n")
                line = line.split("\t")
                self.figures.append(line[0])
                if len(line) > 0:
                    self.captions.append(line[1])
                else:
                    self.captions.append("")

    def modify_captions(self):
        for i in range(len(self.captions)):
            self.captions[i] = re.sub(r"\. ([a-z][-,]*[a-z]*\.)",
                                      r". \\textbf{\1}", self.captions[i])
            self.captions[i] = re.sub(r"^([A-Za-z0-9 ]*\. [A-Za-z0-9\-\(\) ]*\.)",
                                      r"\\textbf{\1}", self.captions[i])
            self.captions[i] = re.sub(r"^([A-Za-z0-9 ]*\.) ",
                                      r"\\textbf{\1} ", self.captions[i])
            self.captions[i] = self.captions[i].replace("<", "\\textless")
            self.captions[i] = self.captions[i].replace(">", "\\textgreater")
            self.captions[i] = re.sub(r"%", r"\%", self.captions[i])
            print(self.captions[i])

    def write_tex_file(self):
        self.outlines = []
        self.outlines = self.add_tex_header(self.outlines)
        for i in range(len(self.figures)):
            self.outlines = self.add_tex_figure(
                self.figures[i], self.captions[i], self.outlines)
        self.outlines.append('\\end{document}')
        with open(self.outtex, 'w') as f:
            for l in self.outlines:
                f.write(l + "\n")

    # TODO: Avoid options hard-coding.
    # TODO: Give option to feed in a template header
    def add_tex_header(self, lines):
        headerlines = ['\\documentclass[11pt]{article}',
                       '\\usepackage[tmargin=.5in, bmargin=.5in,' +
                       'lmargin=.5in, rmargin=.5in]{geometry}',
                       '\\geometry{letterpaper}',
                       '\\usepackage{graphicx}',
                       '\\maxdeadcycles=1000',
                       '\\begin{document}']
        for l in headerlines:
            lines.append(l)
        return(lines)

    def add_tex_figure(self, figure, caption, lines):
        self.suffixes = [f.split(".")[-1] for f in self.figures]
        figlines = ['\\begin{figure}[ht!]',
                    '\\centering',
                    '\\includegraphics[width=\\textwidth,' +
                    'height=0.9\\textheight,keepaspectratio]{' + figure + '}',
                    # Avoid using caption, which makes its own figure number.
                    '\\\\ ' + caption,
                    '\\end{figure}\n']
        for i, l in enumerate(figlines):
            lines.append(l)
            if (i % 10 == 0):
                lines.append('\\clearpage')
        return(lines)

    def compile_tex_file(self):
        os.system("pdflatex " + self.outtex)


if __name__ == "__main__":
    fire.Fire(collate_figures_pdf)
