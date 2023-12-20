from PIL import Image
from sys import argv

image = argv[1]
output = argv[2] + ".png"
print(output)

im = Image.open(image)
im1 = im.crop((0, 0, 700, 360))
im1.save(output, format="PNG")
