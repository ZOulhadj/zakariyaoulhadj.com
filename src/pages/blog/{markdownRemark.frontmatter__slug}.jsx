import React from "react";

import Base from "../../components/base";
import {graphql, Link} from "gatsby";

import "./index.module.scss";
import {GatsbyImage, getImage} from "gatsby-plugin-image";

export default function BlogPost({ data }) {
    const info = data.markdownRemark;
    const post = info.frontmatter;
    const headings = info.headings;
    const bannerImage = getImage(post.banner);

    return (
        <Base>
            <div className="row g-5">
                <div className="col-md-8">

                    <article className="blog-post">
                        <h2 className="blog-post-title mb-1">{post.title}</h2>
                        <p className="blog-post-meta">{post.date}</p>

                        <div dangerouslySetInnerHTML={{ __html: info.html }} />
                    </article>

                    <nav className="blog-pagination" aria-label="Pagination">
                        <a className="btn btn-outline-primary rounded-pill" href="#">Older</a>
                        <a className="btn btn-outline-secondary rounded-pill disabled">Newer</a>
                    </nav>

                </div>

                <div className="col-md-4">
                    <div className="position-sticky" style={{top: 2 +'rem'}}>
                        <div className="p-4">
                            <h4 className="fst-italic">Table of contents</h4>
                            <ol className="list-unstyled mb-0">
                                { headings.map(heading => (
                                    <li className={`ps-0`}>
                                        <Link to={`#${heading.id}`} className="text-dark">{heading.value}</Link>
                                    </li>
                                ))}
                            </ol>
                        </div>

                    </div>
                </div>
            </div>



        {/*    <div className="col-md-8">
                <div className="border rounded p-5 mb-4">
                    <h1 className="fw-bold ps-">{post.title}</h1>
                    <div className="d-flex flex-column text-muted">
                        <span>Posted on {post.date}</span>

                        <span>
                        Time to read: { info.timeToRead } { info.timeToRead > 1 ? "minutes" : "minute" }
                    </span>

                    </div>

                    <GatsbyImage alt={""} image={bannerImage} />
                </div>

                <div dangerouslySetInnerHTML={{ __html: info.html }} />
            </div>

            <div className="col-md-4">
                { headings.map(heading => (
                    <li className={`ps-${heading.depth - 1}`}>
                        <Link to={`#${heading.id}`}>{heading.value}</Link>
                    </li>
                ))}
            </div>*/}
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
      banner {
        childImageSharp {
          gatsbyImageData
        }
      }
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