// Vector saves w, h, d, c, n, k, filter_w(s), filter_h(r), filter_d, pad_w, pad_h, pad_d, wstride, hstride, dstride
std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int,
                       unsigned int, unsigned int, unsigned int, unsigned int,  unsigned int, unsigned int,
                       unsigned int, unsigned int, unsigned int>> training_set = {
    std::make_tuple(2, 2, 2, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(4, 4, 4, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(6, 6, 6, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(8, 8, 8, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(16, 16, 16, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(12, 12, 12, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(24, 24, 24, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(32, 32, 32, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(32, 32, 32, 4, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(48, 48, 48, 4, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(64, 64, 64, 4, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
};

// Vector saves w, h, d, c, n, k, filter_w(s), filter_h(r), filter_d, pad_w, pad_h, pad_d, wstride, hstride, dstride
std::vector<std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int,
                       unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int,
                       unsigned int, unsigned int, unsigned int>> inference_server_set = {
    std::make_tuple(2, 2, 2, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(4, 4, 4, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(6, 6, 6, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(8, 8, 8, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(16, 16, 16, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(12, 12, 12, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(24, 24, 24, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(32, 32, 32, 64, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(32, 32, 32, 4, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(48, 48, 48, 4, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
    std::make_tuple(64, 64, 64, 4, 1, 64, 3, 3, 3, 1, 1, 1, 2, 2, 2),
};

