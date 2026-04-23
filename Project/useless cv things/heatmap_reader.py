import cv2
import numpy as np
import matplotlib.pyplot as plt

# Let's generate a full color synthetic heatmap (Jet colormap simulation)
def create_color_heatmap():
    size = 500
    img = np.zeros((size, size), dtype=np.uint8)
    cv2.circle(img, (150, 150), 80, 200, -1)
    cv2.circle(img, (350, 300), 100, 255, -1)
    cv2.circle(img, (100, 400), 50, 150, -1)
    img = cv2.GaussianBlur(img, (101, 101), 0)
    # Apply JET colormap to make it look like a real color elevation map
    color_img = cv2.applyColorMap(img, cv2.COLORMAP_JET)
    cv2.imwrite('heatmap_color.png', color_img)
    return color_img

#create_color_heatmap()

def extract_color_contours(image_path):
    # 1. Load the COLOR image
    img = cv2.imread(image_path)
    if img is None:
        raise FileNotFoundError(f"Image not found at {image_path}")
    
    # 2. Convert to HSV space (better for color isolation)
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    
    # 3. Define the 'Heat Range' you want to detect.
    # We want to catch from Yellow-Orange to Red (obstacles).
    # NOTE: Red in HSV wraps around 0 and 180, so we use a very broad low range.
    # Hue (color), Saturation (intensity), Value (brightness)
    lower_heat = np.array([0, 150, 50])    # Deep Orange/Red
    upper_heat = np.array([20, 255, 255])  # Pure Yellow
    
    # 4. Create a mask: only pixels within this range will be white.
    binary_mask = cv2.inRange(hsv, lower_heat, upper_heat)
    
    # 5. Extract contours from the mask
    contours, _ = cv2.findContours(binary_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    return contours, img, binary_mask

# --- EXECUTION ---
# You don't pass a single threshold number anymore, the range is defined inside.
curves, color_img, mask = extract_color_contours('heatmap.png')

# --- PLOTTING ---
plt.figure(figsize=(15, 5))

# Plot 1: Source
plt.subplot(1, 3, 1)
plt.imshow(cv2.cvtColor(color_img, cv2.COLOR_BGR2RGB)) # Convert back to RGB for matplotlib
plt.title('Source COLOR Heatmap')
plt.axis('off')

# Plot 2: The Mask (What the drone 'sees' as an obstacle)
plt.subplot(1, 3, 2)
plt.imshow(mask, cmap='gray')
plt.title('Obstacle Mask (Yellow-Red zones)')
plt.axis('off')

# Plot 3: Result
plt.subplot(1, 3, 3)
plt.imshow(cv2.cvtColor(color_img, cv2.COLOR_BGR2RGB))
for i, contour in enumerate(curves):
    pts = contour.reshape(-1, 2)
    pts = np.vstack([pts, pts[0]])
    # Cyan is a good contrast color against JET map
    plt.plot(pts[:, 0], pts[:, 1], color='cyan', linewidth=2.5)

plt.title(f'Detected Borders ({len(curves)} obstacles)')
plt.axis('off')

plt.tight_layout()
plt.show()