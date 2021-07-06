with open('pyne.cpp', 'r', encoding='utf-8') as reader:
    lines = reader.readlines()

writer = open('pyne.cpp.new', 'w', encoding='utf-8')
write_line = True
for i, line in enumerate(lines):
    if 'pyne::Material pyne::Material::decay(' in line:
        write_line = False
        writer.write(line)
        writer.write('  throw pyne::ValueError("Material::decay is not supported in this amalgamated"\n')
        writer.write('                         "version of PyNE.");\n')
    if 'pyne::Material pyne::Material::cram(' in line:
        write_line = False
        writer.write(line)
        writer.write(lines[i + 1])
        writer.write('  throw pyne::ValueError("Material::cram is not supported in this amalgamated"\n')
        writer.write('                         "version of PyNE.");\n')
    if line.strip() == '}':
        write_line = True
    if write_line:
        writer.write(line)
writer.close()
