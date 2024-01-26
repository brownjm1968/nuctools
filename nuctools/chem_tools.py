import re

__all__ = ['get_single_elem_list']

def get_single_elem_list(unary_salts):
    """
    Assuming only chloride salts, break list of unary salt 
    mixtures into a single list of unique chemical elements

    Parameters
    ----------
    unary_salts : list
        List of strings that describe formula for unary salts.

    Returns
    -------
    unique_elements : list
        A list of strings that are the unique elements for all
        salts listed 

    Notes
    -----
    Second author: Tarek Ghaddar

    """
    unique_elements = []
    for salt in unary_salts:
        elem_list = re.findall('[A-Z][^A-Z]*',salt)
        #print(elem_list)
        for elem in elem_list:
            # Getting only the letters from the string
            elem_letters = re.sub('[^a-zA-Z]+', '',elem)
            if elem_letters not in unique_elements:
                unique_elements.append(elem_letters)

    return unique_elements



