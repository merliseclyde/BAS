library(BAS)
library(MASS)
data("UScrime")
crime.bic <- bas.lm(log(y) ~ log(M) + So + log(Ed) +
                             log(Po1) + log(Po2) +
                             log(LF) + log(M.F) + log(Pop) + log(NW) +
                             log(U1) + log(U2) + log(GDP) + log(Ineq) +
                             log(Prob) + log(Time),
                             data = UScrime, n.models = 2^15, prior = "BIC",
                             modelprior = beta.binomial(1, 1),
                            initprobs = "eplogp", pivot = FALSE
                      )
png(file="BAS-image.png")
image(crime.bic)
dev.off()
# install package
imgurl <- system.file("figures/BAS-image.png", package="BAS")
sysfonts::font_add_google("Zilla Slab", "pf", regular.wt = 500)

hexSticker::sticker(imgurl, package="BAS",s_x=1, s_y=1, s_width = 0.6, s_height = 0.7,
        p_x=1, p_y=1.6,  p_size=38, p_color="black", p_family = "pf", 
        p_fontface = "plain",
        filename="inst/logo/logo.png",
        h_color="black", h_fill="white", white_around_sticker = FALSE,
        dpi = 1000 # higher dpi means higher resolution
)

usethis::use_logo("inst/logo/logo.png")
