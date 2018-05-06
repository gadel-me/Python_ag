from PIL import Image, ImageChops

def trim(im):
    """
    Trim the Image by the color of the background.
    Sources:    https://stackoverflow.com/questions/10615901/trim-whitespace-using-pil

    Input:
        > im    opened Image
    """
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)
