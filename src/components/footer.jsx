import React from "react"
import {Link} from "gatsby";
import {useSiteMetadata} from "../hooks/use-site-metadata";

const icons = [
    {
        name: "github",
        url: "https://github.com/ZOulhadj"
    },
    {
        name: "linkedin",
        url: "https://www.linkedin.com/in/zoulhadj/"
    },
    {
        name: "youtube",
        url: "https://www.youtube.com/channel/UCCWqJcNwly8APdxEsY8tZpw"
    },
    {
        name: "instagram",
        url: "https://instagram.com/ZOulhadj"
    },
    {
        name: "stack-overflow",
        url: "https://stackexchange.com/users/8077582/zoulhadj"
    }
]

export default function Footer() {
    const { siteDomain } = useSiteMetadata();

    return (
        <div className="my-4 border-top">
            <footer className="container d-flex flex-wrap justify-content-between align-items-center py-3">
                <div className="col-md-4 d-flex align-items-center">
                    <a href="/" className="mb-3 me-2 mb-md-0 text-muted text-decoration-none lh-1">
                        <svg className="bi" width="30" height="24">
                            <use xlinkHref="#bootstrap"></use>
                        </svg>
                    </a>
                    <span className="mb-3 mb-md-0 text-muted">
                        &copy; { new Date().getFullYear() }
                    </span>
                </div>

                <div className="d-flex align-items-center">
                    { siteDomain }
                </div>

                <ul className="nav col-md-4 justify-content-end list-unstyled d-flex">
                    { icons.map(icon => (
                        <li className="ms-3"><a className="text-muted" href={icon.url} target="_blank" rel="noopener noreferrer">
                            <i className={`bi bi-${icon.name}`}></i>
                        </a></li>
                    ))}
                </ul>
            </footer>
        </div>
    );
}