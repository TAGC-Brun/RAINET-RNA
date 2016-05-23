# ==============================================================
# This class check if two region are overlapped or not
# =============================================================






class OverlapTwoRegion():
    
    # regios_overlap
    # -------------------
    # 
    # Comparison of the possible overlap between two regions
    # Argument:
    #    - Region A and B: two list representing protein regions [start, end] 
    #       
    # Returns:
    #    - a list containing the overlap position or an empty list
    #
    
    @staticmethod
    def region_overlap(region_A, region_B):
        
        item_A = range(region_A[0], region_A[1]+1)
        item_B = range(region_B[0], region_B[1]+1)
        
        overlap = set(item_A) & set(item_B)
        
        return sorted(list(overlap))