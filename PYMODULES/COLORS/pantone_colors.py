__version__ = "2019-05-07"

rgb_colors = (
    0.003922,
    0.309804,
    0.631373,
    0.745098,
    0.188235,
    0.101961,
    0.352941,
    0.356863,
    0.364706,
    0.968627,
    0.572549,
    0.219608,
    0.992157,
    0.949020,
    0.000000,
    0.588235,
    0.364706,
    0.062745,
    0.705882,
    0.709804,
    0.725490,
    0.168627,
    0.639216,
    0.258824,
    0.996094,
    0.996094,
    0.996094,
    0.945098,
    0.454902,
    0.674510,
    0.000000,
    0.576471,
    0.686275,
    0.611765,
    0.239216,
    0.592157,
    0.490196,
    0.643137,
    0.184314,
    0.900000,
    0.400000,
    0.700000,
    0.309804,
    0.133333,
    0.003922,
    0.554688,
    0.777344,
    0.929688,
    0.141176,
    0.125490,
    0.129412,
    0.862745,
    0.796078,
    0.003922,
    1.000000,
    0.960784,
    0.447059,
    0.282353,
    0.415686,
    0.086275,
    0.023529,
    0.352941,
    0.290196,
    0.000000,
    0.396078,
    0.325490,
    0.003922,
    0.682353,
    0.607843,
    0.000000,
    0.113725,
    0.349020,
    0.003922,
    0.611765,
    0.862745,
    0.215686,
    0.243137,
    0.600000,
    0.415686,
    0.125490,
    0.498039,
    0.866667,
    0.000000,
    0.533333,
    0.501961,
    0.019608,
    0.333333,
    0.568627,
    0.000000,
    0.160784,
    0.639216,
    0.133333,
    0.294118,
    0.858824,
    0.431373,
    0.105882,
    0.894118,
    0.552941,
    0.101961,
)

color_names = (
    "blue",
    "red",
    "gray",
    "orange",
    "yellow",
    "tan",
    "silver",
    "green",
    "white",
    "pink",
    "cyan",
    "purple",
    "lime",
    "mauve",
    "ochre",
    "iceblue",
    "black",
    "yellow2",
    "yellow3",
    "green2",
    "green3",
    "cyan2",
    "cyan3",
    "blue2",
    "blue3",
    "violet",
    "violet2",
    "magenta",
    "magenta2",
    "red2",
    "red3",
    "orange2",
    "orange3",
)

pantone_rgb = dict(list(zip(color_names, rgb_colors)))

pantone_colors_hex = {
    "blue": "#3597D8",  # Process Blue
    "red": "#AA2A23",  # 1807
    "gray": "#5B5B5D",  # CG11
    "orange": "#F0903E",  # Orange
    "yellow": "#FAF420",  # Yellow
    "tan": "#915D1B",  # 4635
    "silver": "#B3B3B8",  # CG5
    "green": "#69BE4B",  # 368
    "white": "#FFFFFF",  # standard white
    "pink": "#EB6EAB",  # 224
    "cyan": "#33ADA7",  # 326
    "purple": "#983795",  # Purple
    "lime": "#9CCD67",  # 375
    "mauve": "#D773AE",  # 674
    "ochre": "#79490A",  # 469
    "iceblue": "#8EC7EE",  # 2905
    "black": "#000000",  # standard black
    "yellow2": "#FCF774",  # 3935
    "yellow3": "#BCD943",  # 382
    "green2": "#4B6B1E",  # 371
    "green3": "#A3D49A",  # 359
    "cyan2": "#207066",  # 3292
    "cyan3": "#38B4C6",  # 3125
    "blue2": "#2C72BE",  # 2935
    "blue3": "#255EA8",  # Refl. Blue
    "violet": "#3E3C97",  # Violet
    "violet2": "#681A7E",  # 259
    "magenta": "#E5008B",  # Magenta
    "magenta2": "#E6007C",  # Rubine Red
    "red2": "#E95936",  # warm red
    "red3": "#E7383A",  # Red032
    "orange2": "#D46C27",  # 1525
    "orange3": "#F29855",  # 164
}


def pantone_colors_rgb(color):
    """
    """
    red = int(pantone_colors_hex[color][1:3], 16)
    green = int(pantone_colors_hex[color][3:5], 16)
    blue = int(pantone_colors_hex[color][5:7], 16)
    return (red, green, blue)
