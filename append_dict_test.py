
index = 0
dictionary = {}
for sym in ['a', 'b']:
    dictionary[sym] = {}
    print(dictionary)
    for quantity in ['E', 'V']:
        dictionary[sym][quantity]['value'] = index
        index += 1
        print(dictionary)
