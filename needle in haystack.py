# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""





sequence = input("sequence: ")
flag = input("flag: ")


def test(sequence, flag):
    count = 0
    a = sequence.find(flag, 0)
    if a > -1:
        count += 1
    while a > -1:
        a = sequence.find(flag, (a+1))
        if a > -1:
            count += 1
    return count
        

answer = test(sequence, flag)
print(answer)

'fdas'.find







