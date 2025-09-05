import nuctools.chem_tools as ct

def test_get_protons():
    symbols = ['au','Au',"AU",'b','Ta',"U"]
    gold_Z = [79,79,79,5,73,92]

    for i,symbol in enumerate(symbols):
        assert gold_Z[i] == ct.get_protons(symbol)

def test_get_elem_symbol():
    symbols = ["Au",'B','Ta',"U"]
    protons = [79,5,73,92]

    for i,z in enumerate(protons):
        assert symbols[i] == ct.get_elem_symbol(z)

