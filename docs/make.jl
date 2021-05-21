push!(LOAD_PATH,"../src/")
using Documenter, SignalingToExperts

makedocs(modules = [SignalingToExperts], sitename = "SignalingToExperts.jl")

deploydocs(repo = "github.com/fatimelbah/SignalingToExperts.jl.git", devbranch = "main")
