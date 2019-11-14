import sys
import imageio

s, image = sys.argv

im = imageio.imread(image)
width  = im._meta.EXIF_MAIN.ExifImageWidth
height = im._meta.EXIF_MAIN.ExifImageHeight
xi = 25
x_pixels = [xi]

for _ in range(27):
    xi += 54
    x_pixels.append(xi)

palette_start = 160
yi = 17 + palette_start
y_pixels = [yi]

for i in range(1, 31):

    if i % 6 == 0:
        yi += 41 + 4
    else:
        yi += 41

    y_pixels.append(yi)

pantone_colors = []

for x in x_pixels:
    for y in y_pixels:
        pixel = tuple(im[y][x])
        # skip white and black pixels
        if tuple(pixel) == (255, 255, 255) or tuple(pixel) == (0, 0, 0):
            continue
        else:
            pantone_colors.append(pixel)

        print((x, y, pixel))

#print(pantone_colors)
print((len(set(pantone_colors))))
