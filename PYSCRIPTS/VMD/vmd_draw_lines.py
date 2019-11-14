import graphics
#import molecule

#molecule.load("graphics", "test")
a = 51.2846
b = 46.9652
c = 44.8244
sub_box_length_a = 3.41897333333
sub_box_length_b = 3.35465714286
sub_box_length_c = 3.44803076923
ra = 15
rb = 14
rc = 13
graph_id = 0

for i in range(rb):
    # xy plane
    v1 = (0, i*sub_box_length_b, 0)
    v2 = (a, i*sub_box_length_b, 0)
    graphics.line(graph_id, v1, v2)

    for z in range(rc):
        v1 = (0, i*sub_box_length_b, z*sub_box_length_c)
        v2 = (a, i*sub_box_length_b, z*sub_box_length_c)
        graphics.line(graph_id, v1, v2)

for i in range(ra):
    # xz plane
    v1 = (i*sub_box_length_a, 0, 0)
    v2 = (i*sub_box_length_a, 0, c)
    graphics.line(graph_id, v1, v2)

    for y in range(rb):
        v1 = (i*sub_box_length_a, y*sub_box_length_b, 0)
        v2 = (i*sub_box_length_a, y*sub_box_length_b, c)
        graphics.line(graph_id, v1, v2)


for i in range(rc):
    # yz plane
    v1 = (0, 0, i*sub_box_length_c)
    v2 = (0, b, i*sub_box_length_c)
    graphics.line(graph_id, v1, v2)

    for x in range(ra):
        # yz plane
        v1 = (x*sub_box_length_a, 0, i*sub_box_length_c)
        v2 = (x*sub_box_length_a, b, i*sub_box_length_c)
        graphics.line(graph_id, v1, v2)
