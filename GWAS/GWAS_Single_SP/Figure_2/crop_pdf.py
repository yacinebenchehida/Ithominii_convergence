import sys
import os

def convert_image_to_pdf(image_path, output_path):
    # Check if the input file exists
    if not os.path.exists(image_path):
        print("Error: Input image file not found.")
        return

    # Read input PDF
    with open(image_path, "rb") as image_file:
        reader = PdfFileReader(image_file)
        if reader.numPages == 0:
            print("Error: Input PDF is empty.")
            return
        image = reader.getPage(0)
        # Set lower-left corner of mediabox
        image.mediaBox.lowerLeft = [0, 110]

        # Write to output PDF
        writer = PdfFileWriter()
        writer.addPage(image)
        with open(output_path, "wb") as output_file:
            writer.write(output_file)

if __name__ == "__main__":
    # Check if correct number of arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python script.py input_image output_pdf")
    else:
        input_image = sys.argv[1]
        output_pdf = sys.argv[2]
        convert_image_to_pdf(input_image, output_pdf)
