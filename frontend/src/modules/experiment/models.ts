export interface IExperiment {
    id: number;
    name: string;
    description: string;
    rootDirectory: string;
    location: string;
    meta: object;
    created_at: Date
}

export interface IExperimentUpdate {
    name?: string;
    description?: string;
}

export interface IExperimentCreate {
    name: string;
    rootDirectory: string;
    description?: string;
}
