# PDF
from PyPDF2 import PdfWriter, PdfReader
from sys import argv

image2 = argv[1]
output2 = argv[2] + ".pdf"
reader = PdfReader(image2)
writer = PdfWriter()
image = reader.pages[0]
image.mediabox.lower_left = [0, 110]
writer.add_page(image)

with open(output2, "wb") as fp:
    writer.write(fp)
