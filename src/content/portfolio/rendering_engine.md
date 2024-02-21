---
title: 'Simple Engine (SE)'
description: 'My custom rendering engine built from scratch for applications and experimentation.'
featured: true
category: software
pubDate: 'Sep 10 2023'
updatedDate: 'Nov 22 2023'
heroImage: '@assets/images/portfolio/rendering_engine.png'
---

# About
This is my long-term project that I work on. A general purpose
rendering/simulation engine that I can use to develop various
different applications. Examples include games, simulations,
visualisations etc.

The architecture is design such that the engine and application
are seperate and communicate with each other using a simple
API. The engine is the underlying executable and the application
is a module (dynamic library) which gets loaded during initialisation.
This is so that the engine handles all the platform/renderer specific
functionaility. The benefit of this is that the application will have
no knowledge of specific drivers, platform or renderer APIs and can
instead just contain generic code. This also means that the engine
can implement low-level functionaility without effecting the application
code.

During development, because the application is a dynamic library, it can
be hot-reloaded to quickly see code changes without needing to restart
the entire program.

In terms of memory, the engine will allocate a specific size and a pointer to
this block of memory is passed to the application. The engine itself does not
care how this memory is used. The application uses this block of memory by
casting it to a type that it is aware of. for example "struct application_state;".

This is done so that the game does not need to keep allocating memory all over the place
for instance, when loading a model, or a texture or a large file. Instead, by using a
custom memory allocator, the application can nicely pack memory in which ever way it
wants. Freeing this memory when terminating is even easier, the engine just
calls "free()" once. Thats it.


# Supported platforms
- Windows

# Rendering backends
- OpenGL
- DirectX 11
- Software Raytracing

# Rendering
2D and 3D support
