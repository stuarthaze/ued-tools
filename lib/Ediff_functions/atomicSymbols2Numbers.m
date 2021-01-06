function atomicNumbers = atomicSymbols2Numbers(atomSymbols)
atomicNumbers(strcmp('H',atomSymbols)) = 1;
atomicNumbers(strcmp('C',atomSymbols)) = 6;
atomicNumbers(strcmp('N',atomSymbols)) = 7;
atomicNumbers(strcmp('O',atomSymbols)) = 8;
atomicNumbers(strcmp('P',atomSymbols)) = 15;
atomicNumbers(strcmp('S',atomSymbols)) = 16;
atomicNumbers(strcmp('Fe',atomSymbols)) = 26;
atomicNumbers(strcmp('Se',atomSymbols)) = 34;
atomicNumbers(strcmp('Sn',atomSymbols)) = 50;
atomicNumbers(strcmp('I',atomSymbols)) = 53;
atomicNumbers(strcmp('Pt',atomSymbols)) = 78;
end