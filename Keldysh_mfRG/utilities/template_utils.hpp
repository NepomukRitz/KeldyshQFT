#ifndef FPP_MFRG_TEMPLATE_UTILS_H
#define FPP_MFRG_TEMPLATE_UTILS_H

/// bool pack:
template<bool...> struct bool_pack;
template<bool... bs>
using all_true = std::is_same<bool_pack<bs..., true>, bool_pack<true, bs...>>;

template<class... Ts>
using are_all_integral = all_true<std::is_integral<Ts>::value...>;




template<typename Base, typename T>
constexpr bool is_instance_of(const T*) {
    return std::is_base_of<Base, T>::value;
}


#endif //FPP_MFRG_TEMPLATE_UTILS_H
