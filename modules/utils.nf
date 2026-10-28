def get_species_name(name) {
  def m = name.tr('[A-Z]','[a-z]').tr('.','v').split('/')[0]
  return m
}

def get_gca(name) {
  def m = name.tr('[A-Z]', '[a-z]').tr('.', 'v').replaceAll("_","").split('/').getAt(1)
  return m
}

def concatString(string1, string2, string3){
  return string1 + '_'+string2 + '_'+string3
}






