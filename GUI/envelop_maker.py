from shapely.geometry.polygon import Polygon as Shape
from shapely.ops import unary_union
import shape_sum

def make(l, b):
    x_envelope = []
    y_envelope = []
    
    for i in range(len(l)):
        poly = []
        x_envelope.append([])
        y_envelope.append([])
        
        for k in range(i, len(l)):
            # Safety check: Prevent IndexError out-of-bounds caused by k+1
            # Note: The original code calculated tx, ty but never used them. 
            # The call is preserved here in case it relies on some internal side effects, 
            # though it might be redundant.
            if k + 1 < len(l):
                tx, ty = shape_sum.shape_sum(l[k][:], l[k+1][:], b[k][:], b[k+1][:])
            
            # Use zip to efficiently assemble coordinate pairs
            # Apply .buffer(0) to automatically fix any micro self-intersections in the polygon
            poly.append(Shape(list(zip(l[k], b[k]))).buffer(0))
            
        # Merge all polygons to get the outer envelope
        # Apply .buffer(0) after union to clean up any remaining topological artifacts
        poly_sum = unary_union(poly).buffer(0)
        
        # Safety protection: If the union produces fragmented pieces (MultiPolygon), 
        # extract the main envelope with the largest area to avoid crashes
        if poly_sum.geom_type == 'MultiPolygon':
            poly_sum = max(poly_sum.geoms, key=lambda p: p.area)
            
        # Extract the exterior boundary coordinates of the merged geometry
        ex_x, ex_y = poly_sum.exterior.xy
        x_envelope[i] = list(ex_x)
        y_envelope[i] = list(ex_y)
        
    # Original logic: Append the boundary of the last original body to the end of the envelope list
    if len(l) > 0:
        x_envelope.append(list(l[-1]))
        y_envelope.append(list(b[-1]))
        
    return x_envelope, y_envelope