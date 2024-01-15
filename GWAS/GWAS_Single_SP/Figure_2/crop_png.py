from PIL import Image
from sys import argv

image1 = argv[1]
output1 = argv[2] + ".png"

im = Image.open(image1)
im1 = im.crop((0, 0, 700, 360))
im1.save(output1, format="PNG")
