+++
title = "Why I prefer C instead of C++"
author = ["Zakariya Oulhadj"]
date = 2026-10-10
tags = ["c", "cpp"]
draft = false
+++

**Disclaimer**: This blog post is meant to be a light hearted rant/discussion on
some of C++'s pitfalls and is highly opinionated. I'm sure many will attempt to
find flaws in my arguments so take everything with a grain of salt. End of the
day, you can program in any language you like, just remember to keep things
simple.


## Background {#background}

My journey of learning to program started with C++.

Having been exposed to C++ for over 5 years, my thinking slowly started to
shift. I started to realise that after all these years my knowledge of C++
whilst decent I'd say I always felt that I didn't know everything.

Over the years, through my own research I stumbled across a wide range of
resources

my never ending on [cppreference](https://en.cppreference.com/index.html)

The following people are well known within the game development community and
are considered highly experienced programmers. It is they who have significantly
influenced my way of thinking about code and programming style:

-   [Casey Muratori](https://en.wikipedia.org/wiki/Casey_Muratori) (Handmade Hero)
-   [Johnathan Blow](https://en.wikipedia.org/wiki/Jonathan_Blow) (CEO of Thekla Inc. and creator of Braid and The Witness)
-   Sean Barret (Author of [STB](https://github.com/nothings/stb) libraries)
-   Ryan Fluerer (Author of RAD Debugger)
-   Eskil Steenberg (YouTuber?)
-   Shawn McGrath (Game developer?)

It was only after I discovered Casey Muratori and his Handmade Hero series in
which he develops his own game from scratch that opened my eyes to a different
way of thinking about and writing code.

Ever since I started, I have always been exposed to C++ and object-oriented
programming which has been presented as this holy grail and taught that the
"correct" and clean way to write code was.

Uncle Bob vs Casey Muratori Clean Code Q&amp;A
<https://github.com/unclebob/cmuratori-discussion>

> "Is the pixel an object or is it a group of objects? Is there a container? Do I
> have to ask a factory to get me a color?"

[Source: Getting rid of the OOP mindset](https://youtu.be/GKYCA3UsmrU?t=140)

I'm certaintly not the first and will not be the last when it comes to critising
C++.

When being taught how to program for the first time, students are often
introduced to C first in order to teach them about memory, allocations,
pointers. It's not long, however, until they move onto to learning C++  C

[Why OOP Is A Nightmare](https://youtu.be/C90H3ZueZMM?si=BYGETZNMBqVTmDd9&t=926)


## The Paradox of Choice: The Real Cost of "Flexibility" {#the-paradox-of-choice-the-real-cost-of-flexibility}

C++ programmers often cite the language's vast feature set as its primary
strength, offering the "flexibility" to solve any problem. In reality, this
flexibility is a double-edged sword due to decades of questionable design
decisions and a committee's rigid adherence to backwards compatibility. The
result is a language that has become a "kitchen sink" of abstractions, imposing
a constant, draining **cognitive load** on the developer.

This is the **Paradox of Choice** applied to systems programming: when a language
provides ten different ways to initialize a single variable or pass a piece of
data, it forces the programmer to expend finite mental energy on the how rather
than the what. Instead of focusing on the actual problem, you are trapped in a
cycle of micro-evaluations—weighing "modern" best practices against "legacy"
realities.

Consider a simple analogy: shopping at a hypermarket with fifty aisles of
cereal. You only need a basic meal, yet you are forced to navigate endless
variations, marketing claims, and "modern" vs. "legacy" packaging. By the time
you reach the checkout, you are experiencing Analysis Paralysis. You have
exhausted your mental capital on a trivial choice, leaving you with less focus
for the rest of your day.

C++ is that hypermarket. Every line of code requires a series of
micro-decisions: Should this be a template? Is a virtual destructor necessary
here? Am I using the 'modern' C++20 way or the 'safe' C++11 way? This creates a
culture of bikeshedding and internal debate that has nothing to do with the
machine or the end user. It is a slow accumulation of "choice fatigue" that
makes a codebase feel heavy and unpredictable.

Don't take my word for it. Lets take a look at of all the main ways we can
initialize a variable in C++ (as of writing).

```cpp
int basic = 10;             // 1. C-style copy initialization
int direct(10);             // 2. Functional/Constructor style
int brace{10};              // 3. Uniform initialization (C++11) - Prevents narrowing
int list = {10};            // 4. Copy list initialization (C++11)
int value_init{};           // 5. Value initialization (Zeroes the memory)

auto deduced = 10;          // 6. Type deduction (C++11)
auto braced_deduce{10};     // 7. Direct braced deduction (C++17)

// The "Trap": These look similar but result in different types
auto x = 10;                // x is an int
auto y{10};                 // y is an int
auto z = {10};              // z is a std::initializer_list<int> !!

// C++20 Designated Initializers (Borrowed from C, but with more rules)
struct Point { int x, y; };
Point p = {.x = 1, .y = 2};
```

(Copy Elision C++17) - Look into this

Not sure about you, but personally...? That's just confusing and uneccesarily
complicated for simply just assigning a single value to a variable.

C, by contrast, functions like a small local shop. It provides a limited set of
tools that have remained largely unchanged for decades. There is no "perfect"
abstraction for every niche scenario, but there is also no choice to be made.
You pick up the tool, you apply it to the hardware, and you move on. In C, the
constraints are a pragmatic mercy; they remove the burden of the language
itself, allowing you to focus entirely on the logic of your program.

The purpose of this post is to give my personal experience of the two languages
and why I choose to use C for most projects.

-   Templates
-   Overloading (Unpredictable code flow)
-   Inheritence (links to context)
-   WYSIWYG (No hidden function calls or other side effects)
-   Context (Needing to understand large parts of a code base)


## Overloading {#overloading}

```cpp
int add(int a, int b);       // (A) Integer version
int add(float a, float b);   // (B) Float version
int add(double a, double b); // (C) Double version

add(5, 10);                  // Calls (A): Exact match
add(2.0, 1.0);               // Calls (C): Exact match

// We can't be sure which overload gets called...
add(2.0f, 10);               // Error: Ambigious: (A) or (B) ?
add(5, 5.0);                 // Error: Ambigious: (A) or (C) ?
```

Whilst this will fail to build as a C++ compiler recognises the ambigiouity, the
point is that we as programmers need to understand the types and order to even
know which function may be called. This becomes increasingly tricky in larger
code bases and lots of different function names.

If we compare this to a C version we get the following which makes all function
calls explicit:

```c
int addi(int a, int b);       // (A) Integer version
int addf(float a, float b);   // (B) Float version
int addd(double a, double b); // (C) Double version

addi(10, 20);                 // "i" for integer
addf(10, 20);                 // "f" for float
addd(10, 20);                 // "d" for double
```

A side benfit to having each function be unqiue is that we can easily grep for
the specific one we are looking for which makes searching a codebase
significantly easier.

Take the following example which performs an operation between two vectors
<span class="inline-src language-cpp" data-lang="cpp">`vec_b`</span> and <span class="inline-src language-cpp" data-lang="cpp">`vec_c`</span> and stores the result in <span class="inline-src language-cpp" data-lang="cpp">`vec_a`</span>. A
challenge to the reader: Is this a pairwise multiplication or a dot product?

```cpp
// C++ - Is this mul or dot? We must now check overloaded operator to find
// out...
vec_a = vec_b * vec_c;

// C - Explicit functions clearly defines operation
vector_mul(vec_a, vec_b, vec_c);
vector_dot(vec_a, vec_b, vec_c);
```

Some argue that this is simply 'bad API design' rather than a flaw of C++.
However, C++ is the only language of the two that permits and encourages this
ambiguity. By allowing the programmer to hide complex, potentially heavy math
operations behind a simple \* operator, C++ prioritizes 'math-like' syntax over
engineering clarity. In C, you don't have the option to be ambiguous with
operators. You must name the operation, which forces clarity by default.

A few times I wanted out to out of mere interest to look at the Linux source
code

Whilst I do admit my knowelege of operating systems is very little I did have
this sudden realisation the code I was looking at whilst I did not understand
from a larger context, the specific lines of code were easy to understand

Having said all this, you may come to the conclusion that all this is a skill
issue. This very well may be the case, and I never claimed to be the best
programmer in fact I'd argue the opposite but I rest my case.


## Inheritence (links to context) {#inheritence--links-to-context}


## Code Flow {#code-flow}

When reading a book we are not expected to constantly be jumping


## Templates {#templates}

Out of all C++ language features, templates has to be by far the number one
reason that made be ultimetly switch to C. From not knowing the underlying type
(by definition), significantly longer compile times and larger executables.

There is this infamous [post](http://archive.md/2014.04.28-125041/http://www.boost.org/doc/libs/1_55_0/libs/geometry/doc/html/geometry/design.html) from the [Boost](https://www.boost.org/) C++ library which takes a simple,
easy to reason about <span class="inline-src language-c" data-lang="c">`distance`</span> function that calculates the Euclidean
distance between two points and turns it into a over-engineered templated mess,
enabled by C++.

Before:

```cpp
// Simple, easy to understand code where we can see what the CPU is doing.
double distance(struct point *a, struct point *b) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;

    return sqrt(dx * dx + dy * dy);
}
```

After:

```cpp
// Now we have geometric objects, dispatching, tag dispatching, different
// strategies
template <typename G1, typename G2, typename S>
double distance(G1 const& g1, G2 const& g2, S const& strategy) {
    return dispatch::distance<typename tag<G1>::type,
                              typename tag<G2>::type,
                              G1,
                              G2,
                              S>::apply(g1, g2, strategy);
}
```

The transparency of the original function is sacrificed for a generic
abstraction that provides no immediate clarity. We've traded a five-line
Euclidean formula for a riddle of template parameters. Does <span class="inline-src language-c" data-lang="c">`::apply`</span>
eventually subtract <span class="inline-src language-c" data-lang="c">`x`</span> from <span class="inline-src language-c" data-lang="c">`y`</span>? Probably, but you’ll have to dig
through layers of Boost's "dispatch" logic just to be sure.

Let's take a look at a more extreme example. [EnTT](https://github.com/skypjack/entt/) is a highly popular ECS
(Entity Component System) library. The library itself is great and so I am not
attempting to bash the library but rather use it as a case study in regards to
code readability. There is a specific function <span class="inline-src language-c" data-lang="c">`group`</span> which according to
the library wiki:

> Groups are meant to iterate multiple components at once and to offer a faster
> alternative to multi type views.

The implementation details is not important. If you can understand the code
below then hats off but personally I still struggle to fully understand it. Is
this a **skill issue**? Maybe, but even so, to understand what this code is doing we
need to first understand templates, variadic expansion, compile-time branching,
fold expressions, etc. All of which distract us from what the code is actually
attempting to do which is merely to perform a cache lookup for an existing
component grouping or, if none exists, allocate and initialize a new one.

```cpp
template<typename... Owned, typename... Get, typename... Exclude>
basic_group<owned_t<storage_for_type<Owned>...>,
            get_t<storage_for_type<Get>...>,
            exclude_t<storage_for_type<Exclude>...>>
group(get_t<Get...> = get_t{}, exclude_t<Exclude...> = exclude_t{}) {
    using group_type = basic_group<owned_t<storage_for_type<Owned>...>,
                                   get_t<storage_for_type<Get>...>,
                                   exclude_t<storage_for_type<Exclude>...>>;
    using handler_type = typename group_type::handler;

    if (auto it = groups.find(group_type::group_id()); it != groups.cend()) {
        return {*std::static_pointer_cast<handler_type>(it->second)};
    }

    std::shared_ptr<handler_type> handler{};
    if constexpr (sizeof...(Owned) == 0u) {
        handler = std::allocate_shared<handler_type>(get_allocator(),
                                                     get_allocator(),
                                                     std::forward_as_tuple(assure<std::remove_const_t<Get>>()...),
                                                     std::forward_as_tuple(assure<std::remove_const_t<Exclude>>()...));
    } else {
        handler = std::allocate_shared<handler_type>(get_allocator(),
                                                     std::forward_as_tuple(assure<std::remove_const_t<Owned>>()...,
                                                                           assure<std::remove_const_t<Get>>()...),
                                                     std::forward_as_tuple(assure<std::remove_const_t<Exclude>>()...));
        ENTT_ASSERT(std::all_of(groups.cbegin(), groups.cend(), [](const auto &data) {
            return !(data.second->owned(type_id<Owned>().hash()) || ...);
        }), "Conflicting groups");
    }

    groups.emplace(group_type::group_id(), handler);

    return {*handler};
}
```

I am sure this code has been extensively tested and is highly optimised,
however, from the point of view of a programmer its a nightmare. This has poor
readability, high cognitive load and requires large amount of outside context to
understand.


## Context {#context}

-   WYSIWYG (No hidden function calls or other side effects)
-   Context (Needing to understand large parts of a code base)


## Cognitive Load {#cognitive-load}

With all this being said, I believe this all comes down to the simple principle
of congitive load and clarity. As programmers, we want to focus and solve the
actual problem at hand. Wether thats implementing new functionaility, debugging
existing code etc. We do not want to have to constantly fight the language
itself


## Final Thoughts {#final-thoughts}

Ultimately, this all comes down to

A quote that I think best explains is phomenon:

A quote from Terry Davis, Creator of Temple OS (Full quote [here](https://www.goodreads.com/quotes/10480697-an-idiot-admires-complexity-a-genius-admires-simplicity-a-physicist)):

> An idiot admires complexity, a genius admires simplicity

If you must use C++, I'd suggest you take a look at [Orthodox C++](https://bkaradzic.github.io/posts/orthodoxc++/) which is
according to the blog post:

> Orthodox C++ (sometimes referred as C+) is minimal subset of C++ that improves
> C, but avoids all unnecessary things from so called Modern C++.
