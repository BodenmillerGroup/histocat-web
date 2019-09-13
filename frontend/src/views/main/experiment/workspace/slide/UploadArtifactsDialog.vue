<template>
  <v-dialog v-model="dialog" persistent max-width="600px">
    <template v-slot:activator="{ on }">
      <v-btn v-on="on" small icon color="grey">
        <v-icon small>mdi-cloud-upload</v-icon>
      </v-btn>
    </template>
    <v-card>
      <v-card-title>Upload Artifacts</v-card-title>
      <v-card-text>
        <v-layout column>
          <v-flex>
            <v-file-input
              v-model="file"
              accept="application/zip"
              label="File input"
              display-size="1000"
              clearable
            ></v-file-input>
          </v-flex>
        </v-layout>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn @click="submit" :disabled="!valid">
          Upload
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-dialog>
</template>

<script lang="ts">
import { datasetModule } from "@/modules/datasets";
import { experimentModule } from "@/modules/experiment";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class UploadArtifactsDialog extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly datasetContext = datasetModule.context(this.$store);

  dialog = false;
  file: File | null = null;

  get valid() {
    return this.file !== null;
  }

  reset() {
    this.file = null;
  }

  cancel() {
    this.dialog = false;
  }

  async submit() {
    if (this.valid) {
      const formData = new FormData();
      const file = this.file as File;
      formData.append("file", file, file.name);
      await this.datasetContext.actions.uploadDataset({
        experimentId: this.experimentContext.getters.activeExperimentId,
        data: formData
      });
      this.dialog = false;
    }
  }
}
</script>
