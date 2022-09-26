import React from "react";
import Base from "../../components/base";
import {graphql, Link} from "gatsby";

export default function Index({ data }) {
    return (
        <Base>
            <h1 className={`fw-bold`}>My blog posts</h1>

            { data.allMdx.nodes.map(node => (
                <div className={`border-bottom mb-2 p-2`}>
                    <div className="d-flex justify-content-between align-items-start">
                        <Link to={`/blog/${node.frontmatter.slug}`} className="mb-0 text-decoration-none">
                            <h4>{node.frontmatter.title}</h4>
                        </Link>

                        <div className="d-flex flex-column">
                            <small className="">{node.frontmatter.date}</small>
                            <Link to={`/blog/${node.frontmatter.slug}/#comments`}>
                                <small className="text-muted"> 10 comments</small>
                            </Link>
                        </div>

                    </div>
                </div>
            ))}

            <nav className="d-flex justify-content-center">
                <ul className="pagination">
                    <li className="page-item"><Link className="page-link" to="#">Previous</Link></li>
                    <li className="page-item"><Link className="page-link" to="#">1</Link></li>
                    <li className="page-item"><Link className="page-link" to="#">2</Link></li>
                    <li className="page-item"><Link className="page-link" to="#">Next</Link></li>
                </ul>
            </nav>

        </Base>
    );
}




export const query = graphql`
query {
  allMarkdownRemark(sort: { fields: frontmatter___date, order: DESC}, limit: 20) {
    nodes {
      frontmatter {
        title
        date(formatString: "MMMM D, YYYY")
        slug
        description
      }
      id
    }
  }
}   
`