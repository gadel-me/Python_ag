import cv2


def find_circles(image, show_result=False):
    """
    """
    # read image through command line
    img = cv2.imread(image)
    gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    retval, thresh = cv2.threshold(gray_img, 127, 255, 0)
    img_contours, _ = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

    if show_result is True:
        cv2.drawContours(img, img_contours, -1, (0, 255, 0))

    midpoints = []

    for c in img_contours:
        # compute the center of the contour
        M = cv2.moments(c)
        cX = int(M["m10"] / M["m00"])
        cY = int(M["m01"] / M["m00"])
        midpoints.append([cX, cY])

        if show_result is True:
            cv2.circle(img, (cX, cY), 7, (255, 255, 108), -1)

    if show_result is True:
        cv2.imshow("Image Contours", img)
        cv2.waitKey(0)

    return midpoints
