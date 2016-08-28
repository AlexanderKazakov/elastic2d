#ifndef MACROS_HPP
#define MACROS_HPP

/* *INDENT-OFF* */
#define DO_ONCE(statement) do { statement; } while(0)

#define GET_MACRO_2(_1, _2, MACRO_NAME, ...) MACRO_NAME
#define GET_MACRO_3(_1, _2, _3, MACRO_NAME, ...) MACRO_NAME
#define GET_MACRO_4(_1, _2, _3, _4, MACRO_NAME, ...) MACRO_NAME


/** Useful for Release mode, when some debugging variables are unused */
#define SUPPRESS_WUNUSED(x) DO_ONCE( if (sizeof(x) == 0) throw )

#endif /* MACROS_HPP */
