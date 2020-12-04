<template>
  <span>
    <v-btn @click="trigger" color="primary" elevation="1" small>
      <v-icon small left>mdi-cloud-upload</v-icon>
      {{ label }}
    </v-btn>
    <input :multiple="multiple" class="visually-hidden" type="file" v-on:change="files" ref="fileInput" />
  </span>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { Component, Emit, Prop, Vue } from "vue-property-decorator";

@Component
export default class UploadButton extends Vue {
  readonly projectsContext = projectsModule.context(this.$store);

  @Prop({ type: String, required: true }) label!: string;
  @Prop({ type: Function, required: true }) upload!: (data: FormData) => void;
  @Prop({ default: false }) multiple!: boolean;

  @Emit()
  async files(e): Promise<FileList> {
    const formData = new FormData();
    const file = e.target.files[0];
    formData.append("file", file, file.name);
    e.target.value = "";
    await this.upload(formData);
    return e.target.files;
  }

  trigger() {
    (this.$refs.fileInput as HTMLElement).click();
  }
}
</script>

<style scoped>
.visually-hidden {
  position: absolute !important;
  height: 1px;
  width: 1px;
  overflow: hidden;
  clip: rect(1px, 1px, 1px, 1px);
}
</style>
