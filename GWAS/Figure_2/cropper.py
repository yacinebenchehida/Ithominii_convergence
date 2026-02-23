from PyPDF2 import PdfFileWriter, PdfFileReader
import sys
import os

def crop_pdf(input_pdf, output_pdf):
	# Check if the input file exists
	if not os.path.exists(input_pdf):
		print("Error: Input PDF file not found.")
		return

	# Read input PDF
	with open(input_pdf, "rb") as input_file:
		reader = PdfFileReader(input_file)
		if reader.numPages == 0:
			print("Error: Input PDF is empty.")
			return

		# Create a new PDF writer
		writer = PdfFileWriter()

		# Crop each page in the input PDF and add to writer
		for page_number in range(reader.numPages):
			page = reader.getPage(page_number)
			# Crop bottom of the page (larger value leads to larger cropping)
			page.mediaBox.lowerLeft = [20, 325]
			page.mediaBox.upperRight = [page.mediaBox.getUpperRight_x() - 20, page.mediaBox.getUpperRight_y() - 30]
			writer.addPage(page)

		# Write to output PDF
		with open(output_pdf, "wb") as output_file:
			writer.write(output_file)

if __name__ == "__main__":
	# Check if correct number of arguments are provided
	if len(sys.argv) != 3:
		print("Usage: python script.py input_pdf output_pdf")
	else:
		input_pdf = sys.argv[1]
		output_pdf = sys.argv[2]
		crop_pdf(input_pdf, output_pdf)
