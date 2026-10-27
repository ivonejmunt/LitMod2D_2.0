from shapely.geometry.polygon import Polygon as Shape
import math

def shape_split(vert, l_split, b_split, l, b):
    ###############
    ### Get the coordinates of the split line drawn by the user
    x_drawn = [vert[k][0] for k in range(len(vert))]
    y_drawn = [vert[k][1] for k in range(len(vert))]
    x_st = x_drawn[0] 
    y_st = y_drawn[0]
    x_end = x_drawn[-1]
    y_end = y_drawn[-1]

    ###############
    ### Calculate where the drawn line lands on the original polygon (find the nearest nodes)
    print("start index for split")
    temp_d_st = [math.sqrt((l_split[i] - x_st)**2 + (b_split[i] - y_st)**2) for i in range(len(l_split))]
    split_st_index = temp_d_st.index(min(temp_d_st))
    
    print("end index for split")
    temp_d_end = [math.sqrt((l_split[i] - x_end)**2 + (b_split[i] - y_end)**2) for i in range(len(l_split))]
    split_end_index = temp_d_end.index(min(temp_d_end))
    print('updating path')
    
    # Calculate the x and y of the insertion points
    # Added boundary protection to prevent index out-of-bounds when adding 1
    if split_st_index + 1 < len(l_split):
        l_insert_st = (l_split[split_st_index] + l_split[split_st_index+1]) / 2.0
        b_insert_st = (b_split[split_st_index] + b_split[split_st_index+1]) / 2.0
    
    l_insert_end = (l_split[split_end_index] + l_split[split_end_index-1]) / 2.0
    b_insert_end = (b_split[split_end_index] + b_split[split_end_index-1]) / 2.0

    new_vert = []
    print("adding in the start from the actual body")
    for k in range(0, split_st_index + 1):
        new_vert.append((l_split[k], b_split[k]))
        
    print("adding the drawn points")
    for f in range(len(x_drawn)):
        new_vert.append((x_drawn[f], y_drawn[f]))
        
    print("adding in the end from the actual body") 
    for k in range(split_end_index, len(l_split)):
        new_vert.append((l_split[k], b_split[k]))
    
    ## Assemble the new split polygon and the original polygon
    # Use .buffer(0) to automatically fix any self-intersections caused by freehand drawing
    temp_1_shape = Shape(new_vert).buffer(0)
    
    # Use zip to prevent out-of-bounds errors due to inconsistent array lengths
    temp_2_shape = Shape(list(zip(l_split, b_split))).buffer(0)
    
    # Core geometric calculation: find the difference
    a = temp_2_shape.difference(temp_1_shape) # Get the remaining part after splitting
    err = a.is_valid
    
    # Safety protection: If the cutting operation causes the geometry to become multiple fragments (MultiPolygon),
    # extract the one with the largest area, otherwise .exterior.xy will cause a crash
    if a.geom_type == 'MultiPolygon':
        a = max(a.geoms, key=lambda polygon: polygon.area)
        
    # Final cleanup to ensure strict topological validity
    a = a.buffer(0)
    if a.geom_type == 'MultiPolygon':
        a = max(a.geoms, key=lambda polygon: polygon.area)
        
    print(list(a.exterior.coords))
    xa, ya = a.exterior.xy                  
    
    # Return standard Python lists for compatibility
    return list(xa), list(ya), a, err