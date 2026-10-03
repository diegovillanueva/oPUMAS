# 🔤🔗 M7: translation/dictionary.csv, oPUMAS.jl and the Fortran sources agree (names_check.jl).

@testset "M7 dictionary and twin names in sync" begin
    errors = name_errors()
    foreach(e -> println("❌ ", e), errors)
    @test isempty(errors)
end
