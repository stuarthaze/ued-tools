function MASK = createMaskFromIMG(IMG,cutoff)
MASK = IMG >= cutoff;