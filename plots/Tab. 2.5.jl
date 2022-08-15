function logfct(n, x)
    if iseven(n)
        return x^(floor(n/2))
    else
        return x^(floor(n/2))*log(x) 
    end
end

trueint = [big(1), big(-1), big(0.5), big(-0.25), 1/big(3), -1/big(9)]

# Huybrechs
xh = [
    big"0.005652228205080097135927256306180286802815381484511642501180735040421745840129801",
    big"0.07343037174265227340615889424583036896952475680282060615488198611705797292476631",
    big"0.2849574044625581537145276023871702103641476232120484776096830434953485539434415",
    big"0.6194822640847783814068089433645450272905870223102585701111114420923426831690395",
    big"0.9157580830046983337846091810052298673782425995948626750455409637651145881619066"
]
wh = [
    big"0.02104694579185462911900268287840808491289287018949655448892001847706831638635594",
    big"0.1307055407444466975910762552074437099836611638130339258639980439299088636525566",
    big"0.2897023016713141568415903510284414599387601678790702683025688431810176042568515",
    big"0.350220370120398710285546803898294144066179811444909016579884938718062170377202",
    big"0.2083248416719858061627839069874126010985059866734902347646281556939430453270388"
] 

# without correction
xw = [
    big"0.005652228205080096985364704714164540288439777097029393362217321025654587376822423",
    big"0.07343037174265227291434548988844891067271491723803790122279703941095197694204213",
    big"0.2849574044625581530801816740443905301852311534104494781752696435050821595572108",
    big"0.6194822640847783809759894292692184817366247303019411789029435086195772377551033",
    big"0.9157580830046983336783290669852861774552768452922622721600339157314744085868128"
] 
ww = [
    big"0.02104694579185462879396279961072009964140543466456105177208257643470164849499678",
    big"0.1307055407444466972949547119851833791496176907570239375080691257874139912703965",
    big"0.2897023016713141568802881767439147746790433013130563928516161615535339405867581",
    big"0.3502203701203987106115802151564421924434360098251708467145999192635541789723049",
    big"0.2083248416719858064192140965037395540864975634401877711536322169607962406782587"
]

# corrected
x = [
    big"0.005652228205080097135927256196733224364461736662521571634310808301956553659434744",
    big"0.07343037174265227340615889388832093647469715917834067035684484541883673822337442",
    big"0.2849574044625581537145276019260509101628963210662763809477575315819894593888129",
    big"0.6194822640847783814068089430513733271561626825848934609549153022972420903233235",
    big"0.9157580830046983337846091809279726334110698024540233984987896162377435354971433"
]
w = [
    big"0.02104694579185462911900268264212980706242530912584843774923279348541444162086709",
    big"0.1307055407444466975910762549921867773164715864654795406970237184043603827930743",
    big"0.289702301671314156841590351056571717447074684999630662451019615095549595934383",
    big"0.3502203701203987102855468041352946381567577080211603909627150348339178095411739",
    big"0.2083248416719858061627839071738170600172707113878809681400088381807577699303166"
]

println("Huybrechs:")
for i = 0:5
    println(Float64(abs(
        (sum([wh[j] * logfct(i, xh[j]) for j = 1:5]) - trueint[i+1]) / trueint[i+1]
    )))
end

println("Quadrature without correction function:")
for i = 0:5
    println(Float64(abs(
        (sum([ww[j] * logfct(i, xw[j]) for j = 1:5]) - trueint[i+1]) / trueint[i+1]
    )))
end

println("Quadrature without correction function:")
for i = 0:5
    println(Float64(abs(
        (sum([w[j] * logfct(i, x[j]) for j = 1:5]) - trueint[i+1]) / trueint[i+1]
    )))
end