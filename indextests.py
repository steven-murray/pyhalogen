pos = []
for ix in range(2):
    pos.append([])
    for iy in range(3):
        pos[ix].append([])
        for iz in range(4):
            pos[ix][iy].append(ix+iy+iz)


print pos
