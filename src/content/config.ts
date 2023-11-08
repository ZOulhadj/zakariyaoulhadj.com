import { defineCollection, z } from 'astro:content';

const portfolio = defineCollection({
    // Type-check frontmatter using a schema
    schema: z.object({
	title: z.string(),
	description: z.string(),
	// Transform string to Date object
	pubDate: z.coerce.date(),
	updatedDate: z.coerce.date().optional(),
	heroImage: z.string().optional(),
        github: z.string().optional()
    }),
});

const blog = defineCollection({
    // Type-check frontmatter using a schema
    schema: z.object({
	title: z.string(),
	description: z.string(),
        author: z.string(),
	// Transform string to Date object
	pubDate: z.coerce.date(),
	updatedDate: z.coerce.date().optional(),
	heroImage: z.string().optional(),
        tags: z.array(z.string())
    }),
});



export const collections = { portfolio, blog };
