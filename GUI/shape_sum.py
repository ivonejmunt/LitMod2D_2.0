from shapely.geometry.polygon import Polygon as Shape
from shapely.ops import unary_union  

def shape_sum(x_envelope_1, x_envelope_2, y_envelope_1, y_envelope_2):
    # Build polygon coordinates using list(zip()).
    # Use buffer(0) to repair minor self-intersections in drawn polygons.
    poly_1 = Shape(list(zip(x_envelope_1, y_envelope_1))).buffer(0)
    poly_2 = Shape(list(zip(x_envelope_2, y_envelope_2))).buffer(0)
    
    polygons = [poly_1, poly_2]
    
    # Combine the selected polygons.
    a = unary_union(polygons).buffer(0)
    
    # If union produces disjoint polygons, the existing implementation uses a fallback.
    # The fallback keeps only the largest polygon; the numerical logic is unchanged.
    if a.geom_type == 'MultiPolygon':
        a = max(a.geoms, key=lambda p: p.area)
        
    # Extract the exterior boundary coordinates.
    t_x, t_y = a.exterior.coords.xy
        
    return list(t_x), list(t_y)