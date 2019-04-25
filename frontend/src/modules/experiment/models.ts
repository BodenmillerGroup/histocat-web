export interface IExperiment {
    id: number;
    name: string;
    description: string;
    root_directory: string;
    location: string;
    meta: object;
    created_at: string;
}

export interface IExperimentUpdate {
    name?: string;
    description?: string;
}

export interface IExperimentCreate {
    name: string;
    description?: string;
}
