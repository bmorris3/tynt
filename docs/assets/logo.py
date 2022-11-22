"""
Generate the tynt logo!
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from PIL import Image

np.random.seed(42)

# Save the logo here:
logo_dir = os.path.dirname(__file__)
uncropped_svg_path = os.path.join(logo_dir, 'logo_uncropped.svg')
cropped_svg_path = os.path.join(logo_dir, 'logo.svg')
png_path = os.path.join(logo_dir, 'logo.png')
ico_path = os.path.join(logo_dir, 'logo.ico')

fig, ax = plt.subplots(figsize=(2, 2))

color='darkslategray'

ax.annotate(r'$\rm \alpha$', (0.12, 0.22), fontsize=160, usetex=True, color=color)

n = 3
for i in range(n):
    y = np.linspace(0.35, 0.45, 15)
    x = 5 * (y - 0.4) ** 2 + 0.3 + 0.03 * i
    ax.plot(x, y, color=color, lw=0.8)
ax.add_patch(plt.Circle((0.365, 0.6), radius=0.04, color=color))
ax.set(
    xlim=[0, 1],
    ylim=[0, 1]
)
ax.axis('off')

savefig_kwargs = dict(
    pad_inches=0, transparent=True
)

fig.savefig(uncropped_svg_path, **savefig_kwargs)

# PNG will be at *high* resolution:
fig.savefig(png_path, dpi=800, **savefig_kwargs)

# This is the default matplotlib SVG configuration which can't be easily tweaked:
default_svg_dims = 'width="144pt" height="144pt" viewBox="0 0 144 144"'

# This is a hand-tuned revision to the SVG file that crops the bounds nicely:
custom_svg_dims = 'width="75pt" height="75pt" viewBox="31 20 100 100"'

# Read the uncropped file, replace the bad configuration with the custom one:
with open(uncropped_svg_path, 'r') as svg:
    cropped_svg_source = svg.read().replace(
        default_svg_dims, custom_svg_dims
    )

# Write out the cropped SVG file:
with open(cropped_svg_path, 'w') as cropped_svg:
    cropped_svg.write(cropped_svg_source)

# Delete the uncropped SVG:
os.remove(uncropped_svg_path)

# Convert the PNG into an ICO file:
img = Image.open(png_path)
img.save(ico_path)