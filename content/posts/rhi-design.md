+++
title = "Designing a Graphics Renderer Hardware Interface"
author = ["Zakariya Oulhadj"]
date = 2026-02-10
tags = ["rendering", "renderer-apis", "rhi"]
draft = true
+++

## Introduction {#introduction}

Before we can talk about Renderer Hardware Interfaces we must first understand
the goal of Graphics APIs.

Unfortunatly, if you want your graphics application to run on many different
devices whether thats on PCs, Consoles, due to the vast number of different
graphics APIs some form of abstraction is required. For instance, take a look at
the table below which shows different graphics APIs and on what platforms they
are supported:

| API                  | Windows | Linux | macOS | Android | iOS | Xbox | PS5 |
|----------------------|---------|-------|-------|---------|-----|------|-----|
| ****OpenGL (4.6)**** | ✅      | ✅    | ⚠️     | ✅ (ES) | ⚠️   | ❌   | ❌  |
| ****DirectX 12****   | ✅      | 🛠️     | ❌    | ❌      | ❌  | ✅   | ❌  |
| ****Vulkan (1.4)**** | ✅      | ✅    | 🛠️     | ✅      | 🛠️   | ❌   | ❌  |
| ****Metal (4.0)****  | ❌      | ❌    | ✅    | ❌      | ✅  | ❌   | ❌  |
| ****AGC / GNM****    | ❌      | ❌    | ❌    | ❌      | ❌  | ❌   | ✅  |
| ****WebGPU****       | ✅      | ✅    | ✅    | ✅      | ✅  | ❌   | ❌  |

<div class="NOTES">

-   ✅ ****Native Support**** (First-class citizen)
-   🛠️ ****Layered/Translation**** (Proton, MoltenVK, DXVK)
-   ⚠️ ****Legacy/Deprecated**** (Not maintained)
-   ❌ ****No Support****

</div>

In my opinion, its a mess and in a ideal world we would only need a single API.
Maybe in the future this could hopefully be the case but as of now we need to
deal with this through abstractions.

{{< figure src="/ox-hugo/xkcd-standards.png" >}}

Reference: <https://xkcd.com/927>

So, here is where a Renderer Hardware Interface comes in. It is an abstraction
used in rendering applications (such as rendering engines) where the goal is to
allow us to utilize multiple rendering APIs depending on the platform without
being tied to a specific one. This, however, is not an easy task and can be
incredibly challenging. In fact, in my quest to make my own game engine, this
has been by far the most challenging aspect.

Existing designs range from using C++'s polymorphism feature—where a base
<span class="inline-src language-c" data-lang="c">`RHIDevice`</span> class is inherited by <span class="inline-src language-c" data-lang="c">`VulkanDevice`</span> or
<span class="inline-src language-c" data-lang="c">`OpenGLDevice`</span>—to more modern, data-oriented approaches.

Common patterns include:

-   **Virtual Function Tables (V-Tables)**: The engine calls a generic interface, and
    the CPU looks up the correct backend function at runtime. While intuitive,
    this can introduce a small amount of overhead and often leads to "leaky
    abstractions" where backend-specific logic starts cluttering the base class.

-   **Compile-time Dispatch (Templates/Header Overrides)**: Selecting the backend at
    compile time. This is high-performance but makes it difficult to switch APIs
    without a full rebuild, and it doesn't solve the problem of managing disparate
    state machines.

Instead, I prefer to utilize the concept of a **Deferred Render Command Buffer**.
This moves away from the immediate-mode wrappers typical of many open-source RHI
implementations and toward the more explicit, command-driven flow used by
libraries like [bgfx](https://github.com/bkaradzic/bgfx) as well as in proprietary engines such as EA's [Halcyon](https://media.contentapi.ea.com/content/dam/ea/seed/presentations/wihlidal-halcyonarchitecture-notes.pdf) and
[Our Machinary](https://ruby0x1.github.io/machinery_blog_archive/post/a-modern-rendering-architecture/). This design ensures that the engine never stalls waiting for the
GPU, allowing the front-end to record a frame's worth of intent while the
backend handles the heavy lifting of API-specific execution at a later point in
time.

@@ - Add image of command buffer recording then replaying.


## What should and shouldn't an RHI have? {#what-should-and-shouldn-t-an-rhi-have}

In many RHI designs, there is a tendency for high-level abstractions to 'leak'
into the lower levels. Concepts like meshes and materials are frequent
offenders, but fundamentally, these are gameplay or engine-side constructs. The
GPU is indifferent to what a 'material' represents; it only understands buffers,
pipelines, and state. When an RHI is littered with mesh or entity data, it
creates tight coupling that hinders portability. A clean interface should treat
the GPU as a raw data processor, keeping the 'what' (engine logic) strictly
separated from the 'how' (hardware execution). This is something I've had to
learn the hard way.

An RHI should consist of very few basic concepts such as the ones listed:

-   Buffer
-   Texture
-   Sampler
-   Shader
-   Bind Group
-   Graphics Pipeline
-   Framebuffer

If you find yourself needing to expose "Frames in Flight" to the RHI then I'd
question if this is really needed. A truly generic RHI should hide the messy
fence-and-semaphore dance entirely, providing the illusion of a linear execution
stream even when the hardware is working asynchronously.

The only 'ordering' the engine should care about is data dependency—specifically
through Barriers. If a Compute Shader is writing to a texture that a Fragment
Shader needs to read, the engine should simply insert a barrier command. The RHI
then translates this into API specific primitives. This ensures the engine
focuses on the logic of the data flow, while the RHI handles the
hardware-specific stalls and transitions required to make that flow safe.


## How I Design my RHI Abstraction {#how-i-design-my-rhi-abstraction}

_**Disclaimer**: This design is an evolving project. While the current
implementation is built on an OpenGL 4.6 backend, the architecture is
specifically engineered to support future Vulkan and DirectX 12
implementations._

The system is built on a strict decoupling of responsibilities, split into two
primary layers: a **Renderer Frontend** and a **Renderer Backend**.

The Frontend acts as the producer; its sole responsibility is to record a
sequence of high-level, API-agnostic commands into a linear array. At the end of
the frame, the Backend takes over as the consumer. It "replays" this command
stream, dispatching each instruction to its API-specific implementation. This
"Record-and-Playback" model ensures that the engine's core logic remains
entirely isolated from the underlying graphics hardware.

{{< figure src="/ox-hugo/rhi-design.png" >}}

At the heart of the RHI is the <span class="inline-src language-c" data-lang="c">`render_commands`</span> structure, which manages
the lifecycle of the frame through two dedicated byte-streams: one for Resources
and one for Rendering.

Rather than using a fixed-size array or a union of commands, I've implemented
these buffers as a contiguous block of bytes. The reasoning here is twofold:

-   **Memory Efficiency**: Commands vary significantly in size. A
    <span class="inline-src language-c" data-lang="c">`cmd_buffer_destroy`</span> might only need a few bytes for a handle, whereas a
    <span class="inline-src language-c" data-lang="c">`cmd_draw`</span> might exceed 64 bytes to accommodate pipelines, buffers, and
    bindings. Using a union would force every command to occupy the space of the
    largest possible type—leading to massive memory fragmentation. By packing
    commands tightly into a byte array, we minimize the memory footprint.

-   **Cache Locality**: A linear byte-stream is incredibly friendly to the CPU cache.
    During the backend's "playback" phase, the CPU can pre-fetch the command data
    sequentially, avoiding the pointer-chasing and cache misses associated with
    linked lists or arrays of pointers.

<!--listend-->

```c
struct command_buffer {
    u8 *buffer;
    u64 capacity;
    u64 offset;

    u32 command_count;
    u64 read_offset;
};

struct render_commands {
    struct command_buffer resource_cmd_buffer;
    struct command_buffer render_cmd_buffer;

    // Pointer to memory mapped GPU staging buffer
    u8 *g_staging_buffer;
    u64 g_staging_size;
    u64 max_staging_size;
};
```


### Commands {#commands}

To facilitate a stateless and deferred architecture, I have defined various
command types. These are categorized into Resource Commands (handling the
lifecycle and data of GPU objects) and Rendering Commands (handling the actual
submission of work to the hardware).

```c
enum render_cmd_type {
    /* Resource commands */
    Command_Buffer_Create,
    Command_Buffer_Destroy,
    Command_Buffer_Upload,
    Command_Buffer_Copy,
    Command_Texture_Create,
    Command_Texture_Destroy,
    Command_Texture_Upload,
    Command_Texture_Copy,
    Command_Sampler_Create,
    Command_Sampler_Destroy,
    Command_Shader_Create,
    Command_Shader_Destroy,
    Command_Bind_Group_Create,
    Command_Bind_Group_Destroy,
    Command_Graphics_Pipeline_Create,
    Command_Graphics_Pipeline_Destroy,
    Command_Framebuffer_Create,
    Command_Framebuffer_Destroy,

    /* Rendering commands */
    Command_Framebuffer_Begin,
    Command_Framebuffer_End,
    Command_Set_Viewport,
    Command_Draw,
    Command_Present,
    Command_Debug_Marker_Begin,
    Command_Debug_Marker_End,
    Command_Timer_Query_Begin,
    Command_Timer_Query_End,
};
```

Every entry in the byte-stream follows a strict protocol consisting of a Header
followed by a Data Payload. The header serves as the "Instruction Pointer" for
the backend, providing the necessary metadata to decode the variably-sized data
that follows.

```c
struct render_cmd_header {
    enum render_cmd_type type; // What is this command?
    u32 size;                  // How many bytes is the total payload?
};
```

By embedding the size in the header, the backend dispatcher knows exactly how
many bytes to skip to reach the next command, allowing the stream to be
traversed linearly.

For instance, a buffer creation command encapsulates the header, the
pre-allocated handle, and the configuration descriptors into a single,
tightly-packed structure:

```c
struct render_cmd_buffer_create {
    struct render_cmd_header header;

    struct buffer_handle buffer;     // The ID reserved by the frontend
    struct buffer_create_info info;  // Usage flags, size, and memory type
};
```

Once the frontend has finished recording the frame's intent, the command buffer
is handed off to the backend for execution. A major advantage of this
architecture is that each backend (OpenGL, Vulkan, Metal) has total autonomy
over how it interprets the stream. If a specific platform doesn't support a
feature—such as certain debug markers or specific timer queries—the backend can
simply ignore those commands without affecting the frontend logic.

In my current OpenGL implementation, the backend processes the buffer by
iterating through the byte-stream and "dispatching" each entry. For every
command encountered, the backend identifies the type via the header and calls
the corresponding API-specific function to update the GPU state or issue a draw
call.

```c
#define COMMAND_ENTRY(type, command, function) \
    case type: {                               \
    const command *cmd = (command *)header;    \
        function(gl, cmd);                     \
    } break                                    \

static void
opengl_handle_commands(struct opengl_backend *gl, struct command_buffer *cb)
{
    while (command_buffer_has_commands(cb)) {
        struct render_cmd_header *header = command_buffer_get_next_command(cb);
        switch (header->type) {
            /* Resource Dispatch */
            COMMAND_ENTRY(Command_Buffer_Create,             struct render_cmd_buffer_create,             opengl_buffer_create);
            COMMAND_ENTRY(Command_Buffer_Destroy,            struct render_cmd_buffer_destroy,            opengl_buffer_destroy);
            COMMAND_ENTRY(Command_Buffer_Upload,             struct render_cmd_buffer_upload,             opengl_buffer_upload);
            COMMAND_ENTRY(Command_Buffer_Copy,               struct render_cmd_buffer_copy,               opengl_buffer_copy);
            COMMAND_ENTRY(Command_Texture_Create,            struct render_cmd_texture_create,            opengl_texture_create);
            COMMAND_ENTRY(Command_Texture_Destroy,           struct render_cmd_texture_destroy,           opengl_texture_destroy);
            COMMAND_ENTRY(Command_Texture_Upload,            struct render_cmd_texture_upload,            opengl_texture_upload);
            COMMAND_ENTRY(Command_Texture_Copy,              struct render_cmd_texture_copy,              opengl_texture_copy);
            COMMAND_ENTRY(Command_Sampler_Create,            struct render_cmd_sampler_create,            opengl_sampler_create);
            COMMAND_ENTRY(Command_Sampler_Destroy,           struct render_cmd_sampler_destroy,           opengl_sampler_destroy);
            COMMAND_ENTRY(Command_Framebuffer_Create,        struct render_cmd_framebuffer_create,        opengl_framebuffer_create);
            COMMAND_ENTRY(Command_Framebuffer_Destroy,       struct render_cmd_framebuffer_destroy,       opengl_framebuffer_destroy);
            COMMAND_ENTRY(Command_Shader_Create,             struct render_cmd_shader_create,             opengl_shader_create);
            COMMAND_ENTRY(Command_Shader_Destroy,            struct render_cmd_shader_destroy,            opengl_shader_destroy);
            COMMAND_ENTRY(Command_Bind_Group_Create,         struct render_cmd_bind_group_create,         opengl_bind_group_create);
            COMMAND_ENTRY(Command_Bind_Group_Destroy,        struct render_cmd_bind_group_destroy,        opengl_bind_group_destroy);
            COMMAND_ENTRY(Command_Graphics_Pipeline_Create,  struct render_cmd_graphics_pipeline_create,  opengl_pipeline_create);
            COMMAND_ENTRY(Command_Graphics_Pipeline_Destroy, struct render_cmd_graphics_pipeline_destroy, opengl_pipeline_destroy);
            COMMAND_ENTRY(Command_Compute_Pipeline_Create,   struct render_cmd_compute_pipeline_create,   opengl_compute_pipeline_create);
            COMMAND_ENTRY(Command_Compute_Pipeline_Destroy,  struct render_cmd_compute_pipeline_destroy,  opengl_compute_pipeline_destroy);

            /* Render Dispatch */
            COMMAND_ENTRY(Command_Framebuffer_Begin,  struct render_cmd_framebuffer_begin,  opengl_framebuffer_begin);
            COMMAND_ENTRY(Command_Framebuffer_End,    struct render_cmd_framebuffer_end,    opengl_framebuffer_end);
            COMMAND_ENTRY(Command_Set_Viewport,       struct render_cmd_viewport,           opengl_set_viewport);
            COMMAND_ENTRY(Command_Draw,               struct render_cmd_draw,               opengl_draw);
            COMMAND_ENTRY(Command_Dispatch,           struct render_cmd_dispatch,           opengl_dispatch);
            COMMAND_ENTRY(Command_Present,            struct render_cmd_present,            opengl_present);
            COMMAND_ENTRY(Command_Debug_Marker_Begin, struct render_cmd_debug_marker_begin, opengl_debug_marker_begin);
            COMMAND_ENTRY(Command_Debug_Marker_End,   struct render_cmd_debug_marker_end,   opengl_debug_marker_end);
            COMMAND_ENTRY(Command_Timer_Query_Begin,  struct render_cmd_timer_query_begin,  opengl_timer_query_begin);
            COMMAND_ENTRY(Command_Timer_Query_End,    struct render_cmd_timer_query_end,    opengl_timer_query_end);
        default: {
            assert(false && "Unknown command type encountered in backend dispatcher");
        } break;
        }
    }
}
```


### Resources Commands {#resources-commands}

My approach to RHI is to take this a step further and have the creation and
deletion of GPU resources be part of the render command system.

Here is a snippet of some functions that are used to push resource commands to
the resource command buffer. Note that <span class="inline-src language-c" data-lang="c">`cb`</span> would just be a pointer to the
<span class="inline-src language-c" data-lang="c">`render_commands`</span> struture but has been omitted for the purposes of
clarity.

```c
void cmd_buffer_create(cb, struct buffer_handle id, const struct buffer_create_info info);
void cmd_buffer_destroy(cb, struct buffer_handle id);
void cmd_buffer_upload(cb, struct buffer_handle id, const void *data, u32 size, u32 dst_offset);
void cmd_buffer_copy(cb, struct buffer_handle src_id, struct buffer_handle dst_id, u32 src_offset, u32 dst_offset, u32 size);
void cmd_texture_create(cb, struct texture_handle id, const struct texture_create_info info);
void cmd_texture_destroy(cb, struct texture_handle id);
void cmd_texture_upload(cb, struct texture_handle id,
                        u32 mip_level,
                        u32 x_offset, u32 y_offset, u32 z_offset,
                        u32 width, u32 height, u32 depth,
                        const void *data, u32 size);
void cmd_texture_copy(cb, struct texture_handle src_id, struct texture_handle dst_id);

// Simmilar commands for samplers, shaders, bind groups, pipelines,
// framebuffers, etc.
```

This raises some interesting and important questions.

-   (1) **Identity: How are handles (IDs) generated and managed?**
-   (2) **Availability: How can the engine reference a GPU resource before the hardware has actually allocated it?**
-   (3) **Ordering: How do we ensure a draw call doesn't attempt to use a resource
    that hasn't finished initializing?**

To solve these, we must shift our mental model of resource management. In this
system, handles are a front-end responsibility. When the engine requires a new
vertex buffer, it doesn't wait for a GPU pointer; it simply reserves a unique ID
from a handle manager. This ID is then passed to the various <span class="inline-src language-c" data-lang="c">`cmd_*`</span>
functions, acting as a "contract" that the backend will eventually fulfill by
associating that ID with a concrete GPU object.

Regarding availability, it is vital to remember that all commands are deferred.
The front-end is not interacting with live GPU memory; it is recording a
manifest of future work. The "resource" exists as a logical concept in the
command stream long before it exists as a physical allocation in VRAM.

Finally, to address ordering, the <span class="inline-src language-c" data-lang="c">`render_commands`</span> structure maintains two
distinct streams: a Resource Command Buffer and a Render Command Buffer. By
segregating "structural" commands (creation, destruction, and uploads) from
"operational" commands (draw calls and state changes), we can guarantee
execution order. The backend always flushes the resource buffer first, ensuring
that every buffer, texture, and pipeline state is fully baked and resident in
memory before the first draw call is ever dispatched.

Let's take a look at simple example for creating a buffer. Below is all the
definitions used.

```c
enum buffer_usage_flags {
    Buffer_Usage_Vertex         = 1 << 0,
    Buffer_Usage_Index          = 1 << 1,
    Buffer_Usage_Uniform        = 1 << 2,
    // Other types can be added ...
};

enum buffer_memory_type {
    Buffer_Memory_Device_Local,
    Buffer_Memory_Host_Visible,
};

struct buffer_create_info {
    enum buffer_usage_flags usage_flags;
    enum buffer_memory_type memory_type;
    u32 size;
};

struct render_cmd_buffer_create {
    struct render_cmd_header header;

    struct buffer_handle buffer;
    struct buffer_create_info info;
};

void cmd_buffer_create(struct render_commands *cb,
                       struct buffer_handle buffer,   // uint16_t ID handle
                       const struct buffer_create_info info);
```

```c
struct buffer_create_info info = {
    .usage_flags = Buffer_Usage_Vertex |
                   Buffer_Usage_Index;
    .memory_type = Buffer_Memory_Device_Local,
    .size = 100
};

struct buffer_handle my_handle = { ... }; // Engine generated uint16_t ID

cmd_buffer_create(cb, my_handle, info);
```

Notice that the buffer creation command intentionally excludes data
initialization. Much like the standard library separates memory allocation
<span class="inline-src language-c" data-lang="c">`malloc`</span> from data population <span class="inline-src language-c" data-lang="c">`memcpy`</span>, the RHI treats resource
residency and data transfer as distinct operations. This is handled via a
dedicated <span class="inline-src language-c" data-lang="c">`cmd_buffer_upload`</span> command.

Because these commands are deferred, we face a synchronization challenge: the
source data must remain valid until the backend eventually processes the queue.
We solve this by copying the user's data into a staging area—either a transient
CPU ring buffer or, for maximum performance, a persistently mapped GPU buffer
that the backend can read from directly during the dispatch phase.


### Render Commands {#render-commands}

Unlike resource commands, there are significantly fewer render-specific commands
as shown below:

```c
void cmd_framebuffer_begin(cb, struct framebuffer_handle fb_id, const f32 clear_color[4], f32 clear_depth, u32 clear_stencil);
void cmd_framebuffer_end(cb, struct framebuffer_handle fb_id);
void cmd_set_viewport(cb, u32 x, u32 y, u32 width, u32 height);
void cmd_draw(cb,
              struct graphics_pipeline_handle pipeline_id,
              const struct bind_group_handle *groups,
              u32 group_count,
              struct buffer_handle *vbs_ids,
              u32 vb_count,
              struct buffer_handle ib_id,
              enum draw_index_type index_type,
              u32 element_count,
              u32 vertex_offset,
              u32 index_offset,
              u32 instance_offset,
              u32 instance_count,
              struct scissor_rect scissor);
```

You'll notice that the <span class="inline-src language-c" data-lang="c">`cmd_draw`</span> function is quite parameter-heavy. This
is a deliberate design choice to move away from stateful rendering in favor of a
stateless architecture. Rather than relying on separate commands to bind buffers
or set pipeline states—which often leads to the 'state leakage' notorious in
OpenGL—every piece of data required to fulfill a draw call is encapsulated
within a single, self-contained command. This ensures that each draw call is
independent; it neither inherits state from the previous command nor leaks its
own state into the next, providing a much higher degree of predictability and
thread-safety.

While the frontend interface is stateless, the backend is responsible for
redundant state filtering. To avoid the heavy performance cost of unnecessary
driver calls, the backend tracks the currently bound resource handles. When a
new draw command is dispatched, the backend compares the requested handles
(pipelines, buffers, bind groups) against the previous ones; if the IDs match,
the bind call is skipped entirely.


### Benefits of this Design {#benefits-of-this-design}

-   **Backend Isolation** - Within an engine's renderer, render backend commands
    issued simply append (writes memory) to a linear buffer with no API-specific
    calls begin made which subsequently means no GPU driver communication. All GPU
    communication only happens when the command stream is being replayed leading
    to predictable communication.

    In fact, as I do within my own engine, since no API-specific calls are made,
    the render backend can be compiled as a completly separate shared library that
    gets loaded at runtime allowing for it to be hot-reloaded or changed
    dynamically to a different backend.

-   **Debugging/Replaying** - Since the entire is render command buffer is recorded,
    it can easily be saved/restored for debugging purposes.

-   **Stateless Rendering** -


## Final Thoughts {#final-thoughts}

I hope this brief dive into my RHI architecture provides some value to others
navigating the same "abstraction hell." While the system is still very much a
work in progress, and certain specifics regarding explicit APIs like Vulkan or
DX12 are still being refined, I am confident that this Deferred Command approach
is the right foundation. It moves the complexity away from the engine logic and
puts it exactly where it belongs: in the backend translator.

Designing a renderer is rarely about finding the "perfect" solution on the first
try; it is about building a contract between your engine and the hardware that
is flexible enough to evolve. As I continue to battle with the nuances of modern
graphics programming, I will post further updates and deep dives into major
architectural changes.

If you're building your own engine, my advice is simple: keep your RHI "dumb,"
your commands stateless, and your data local. The GPU doesn't care about your
game's abstractions—it just wants its buffers and its state.


## Helpful Resources {#helpful-resources}

-   <https://www.gijskaerts.com/wordpress/?p=112>
-   <https://media.contentapi.ea.com/content/dam/ea/seed/presentations/wihlidal-halcyonarchitecture-notes.pdf>
-   <https://alextardif.com/RenderingAbstractionLayers.html>
-   <https://www.sebastianaaltonen.com/blog/no-graphics-api>
-   <https://www.gamedev.net/forums/topic/694968-renderer-design-bitquid-our-machinery/>
-   <https://zeux.io/2020/02/27/writing-an-efficient-vulkan-renderer/>
-   <https://advances.realtimerendering.com/s2023/AaltonenHypeHypeAdvances2023.pdf>
-   <https://bitsquid.blogspot.com/2017/02/stingray-renderer-walkthrough.html>
-   <http://alinloghin.com/articles/command_buffer.html>
-   <https://blog.mecheye.net/2023/09/how-to-write-a-renderer-for-modern-apis/>
-   <https://ruby0x1.github.io/machinery_blog_archive/post/a-modern-rendering-architecture/>
