import { defineCollection, z } from "astro:content";

const portfolioCollection = defineCollection({
    type: "content",
    schema: ({ image }) => z.object({
        title: z.string(),
        description: z.string(),
        featured: z.boolean().default(false),
        // Transform string to Date object
        pubDate: z.coerce.date(),
        updatedDate: z.coerce.date().optional(),
        heroImage: image(),
        website: z.string().optional(),
        github: z.string().optional()
    })
});

const blogCollection = defineCollection({
    type: "content",
    schema: ({ image }) => z.object({
        title: z.string(),
        description: z.string(),
        author: z.string(),
        // Transform string to Date object
        pubDate: z.coerce.date(),
        updatedDate: z.coerce.date().optional(),
        heroImage: image(),
        draft: z.boolean()
    })
});

export const collections = {
    "portfolio": portfolioCollection,
    "blog": blogCollection,
};
