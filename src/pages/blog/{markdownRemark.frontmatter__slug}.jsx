import React from "react";

import Base from "../../components/base";
import {graphql, Link} from "gatsby";

import "./index.module.scss";

export default function BlogPost({ data }) {
    const info = data.markdownRemark;
    const post = info.frontmatter;
    const headings = info.headings;

    return (
        <Base>
            <div className="mb-4">
                <h1 className="fw-bold ps-">{post.title}</h1>
                <div className="d-flex flex-column text-muted">
                    <span>Posted on {post.date}</span>

                    <span>
                        Time to read: { info.timeToRead } { info.timeToRead > 1 ? "minutes" : "minute" }
                    </span>

                    { headings.map(heading => (
                        <li className={`ps-${heading.depth - 1}`}>
                            <Link to={`#${heading.id}`}>{heading.value}</Link>
                        </li>
                    ))}
                </div>

            </div>

            <div dangerouslySetInnerHTML={{ __html: info.html }} />
        </Base>
    );
}

export const query = graphql`
query ($id: String) {
  markdownRemark(id: {eq: $id}) {
    frontmatter {
      title
      date(formatString: "MMMM D, YYYY")
      description
    }
    headings {
      id
      value
      depth
    }
    html
    timeToRead
  }
}    
    
`