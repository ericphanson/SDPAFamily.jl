function install_custom_WSL_binaries(prefix, verbose)
    prefix = Prefix(joinpath(prefix, "bin"))
    BinaryProvider.download_verify("https://github.com/ericphanson/SDPA_GMP_Builder/raw/master/deps/sdpa_gmp_wsl",  # url
            "f8ed0c3f2aefa1ab5a90f1999c78548625f6122f969972b8a51b54a0017b3a59", # hash
            joinpath(prefix.path,"sdpa_gmp");  # destination
            verbose=verbose)

    products = Product[
        FileProduct(prefix, "sdpa_gmp", :sdpa_gmp),
    ]

    return products

end
