from shapely.geometry.polygon import Polygon as Shape

def shape_substract(x_envelope_1, x_envelope_2, y_envelope_1, y_envelope_2):
    # Use list(zip()) to efficiently assemble coordinate pairs
    # Add .buffer(0) to automatically fix any micro self-intersections or drawing artifacts
    poly_1 = Shape(list(zip(x_envelope_1, y_envelope_1))).buffer(0)
    poly_2 = Shape(list(zip(x_envelope_2, y_envelope_2))).buffer(0)
    
    # Geometric calculation: subtract poly_1 from poly_2
    a = poly_2.difference(poly_1)
    
    # Final cleanup to ensure strict topological validity after the difference operation
    a = a.buffer(0)
    
    # Safety protection: If the subtraction causes the polygon to break into multiple fragments (MultiPolygon)
    # extract the one with the largest area to prevent crashes when calling .exterior later
    if a.geom_type == 'MultiPolygon':
        a = max(a.geoms, key=lambda polygon: polygon.area)
        
    t_x, t_y = a.exterior.xy
    
    # Return standard Python lists for compatibility
    return list(t_x), list(t_y)