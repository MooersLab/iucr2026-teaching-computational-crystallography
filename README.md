![Version](https://img.shields.io/static/v1?label=sofware-fayre-2026-slides&message=0.1&color=brightcolor)
![License](https://img.shields.io/badge/License-MIT-yellow)
![License](https://img.shields.io/badge/License-CC4-orange)

## Under construction; COME BACK Later!

## What is this?
Slideshow for the XXVII Congress of the International Union of Crystallography (IUCr) 2026 talk in the microsymposium MS-116, Teaching Computational Crystallography, on August 17, 2026, from 3:00 to 3:30 PM in Room 224 of the BMO Convention Center in Calgary, Alberta. 
This repo includes the Jupyter Notebooks used in the talk.
The repo also includes the source files for assembly of the PDF with LaTeX.
Run `make` in the top folder, provided TeX Live or MiKTeX is already installed.

## Compile the slides yourself

The `main.tex` file is the source file for the slides.
The `main.pdf` contains the slides. 
The `main.tex` is used with LaTeX and Beamer to generate the slides.
The `Makefile` generates the PDF when `make` is run.
It opens with the resulting PDF in Skim.app when you enter `make view`.

## LICENSE

The source code is covered by the MIT License.
The images are covered by a Creative Commons Version 4 license.


## Funding

- NIH: R01 CA242845.
- NIH: P30 CA225520 (PI: R. Mannel).
- NIH: P20 GM103640 and P30 GM145423 (PI: A. West).
