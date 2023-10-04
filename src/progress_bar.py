# progress bar code credit: https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters

def progress_bar(iteration, total, step, nrg, atoms, prefix = '', suffix = '', length = 100, fill = '█', up = '\x1B[3A', clr = '\x1B[0K'):
   percent = ('{0:.' + '1' + 'f}').format(100 * (iteration / float(total)))
   filled_length = int(length * iteration // total)
   bar = fill * filled_length + '-' * (length - filled_length)
   print(f'{up}{prefix} |{bar}| {percent}% {suffix}{clr}\nStep: {step}, Potential Energy: {nrg} eV, Atoms: {atoms}{clr}\n')

# print new line when complete

   if iteration == total: 
      print('\nCOMPLETE!')