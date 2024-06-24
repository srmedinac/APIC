from histoprep import SlideReader

# Read slide image.
reader = SlideReader("/path/to/slide.svs")
# Detect tissue.
threshold, tissue_mask = reader.get_tissue_mask(level=-1)
tile_coordinates = reader.get_tile_coordinates(
    tissue_mask, width=1024, overlap=0, max_background=0.7
)
# Save patches
tile_metadata = reader.save_regions(
    "/save/path", tile_coordinates, threshold=threshold, save_metrics=True
)
