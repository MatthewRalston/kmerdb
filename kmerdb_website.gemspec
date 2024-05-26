# coding: utf-8

Gem::Specification.new do |spec|
  spec.name          = "kmerdb_jekyll_theme_forty"
  spec.version       = "1.2"
  spec.authors       = ["Matthew Ralston"]
  spec.email         = ["mralston.development@gmail.com"]

  spec.summary       = %q{A website for bioinformatics computation based on the forty theme by HTML5 UP.}
  spec.homepage      = "https://github.com/MatthewRalston/kmerdb"
  spec.license       = "Apache"

  spec.files         = `git ls-files -z`.split("\x0").select { |f| f.match(%r{^(assets|_layouts|_includes|_sass|LICENSE|README)}i) }

  spec.add_development_dependency "jekyll", "~> 4.2"
  spec.add_development_dependency "bundler", "~> 2.5"
  spec.add_development_dependency "webrick", "~> 1.8"
end
