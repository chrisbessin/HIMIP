import os, sys

from docx import Document
from docx.shared import Inches


document = Document()

document.add_heading('Tests methods samples', 0)

dir = 'D:\\SOLSA\\HIMIP\\20191007_solsa_himip\\interface'

os.listdir(dir)

document.save('demo.docx')
